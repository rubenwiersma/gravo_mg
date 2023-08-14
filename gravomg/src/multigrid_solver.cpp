#define _USE_MATH_DEFINES
#include "gravomg/multigrid_solver.h"
#include "gravomg/utility.h"
#include "gravomg/sampling.h"

#include <igl/invert_diag.h>

#include <cmath>
#include <numeric>
#include <chrono>

#include <Eigen/Eigenvalues>

namespace MGBS {
	/* Constructor */
	MultigridSolver::MultigridSolver(
		Eigen::MatrixXd& V, Eigen::MatrixXi& neigh, Eigen::SparseMatrix<double>& M) : 
		V(V), neigh(neigh), M(M) {
		igl::invert_diag(M, Minv);
		V0 = V;
		hierarchyTiming["n_vertices"] = V.rows();
	}

	/* ====================== BASIC FUNCTIONS  ======================= */
	MultigridSolver::~MultigridSolver()
	{
		// std::cout << "Deleting data!\n";
		V.resize(0, 0);
		V0.resize(0, 0);
		neigh.resize(0, 0);
		DoF.clear(); DoF.shrink_to_fit();
		U.clear(); U.shrink_to_fit();
		Abar.clear(); Abar.shrink_to_fit();
		Mbar.clear(); Mbar.shrink_to_fit();
		neighHierarchy.clear(); neighHierarchy.shrink_to_fit();
		nearestSource.clear(); nearestSource.shrink_to_fit();
		geodDistance.clear(); geodDistance.shrink_to_fit();
		samples.clear(); samples.shrink_to_fit();
		PointsToSampleMap.clear(); PointsToSampleMap.shrink_to_fit();
	}

	// Build the hierarchy
	void MultigridSolver::buildHierarchy()
	{
		nearestSource.clear();
		nearestSource.shrink_to_fit();

		plf::nanotimer timer;
		timer.start();
		if (sig06) {
			constructProlongationSIG06();
		} else if (ablation) {
			constructProlongationAblation();
		} else {
			constructProlongation();
		}
		hierarchyTiming["hierarchy"] = timer.get_elapsed_ms();

		UOurs = U;
	}

	void MultigridSolver::constructProlongation()
	{
		// Create prolongation operator for level k+1 to k
		// Points in current level
		Eigen::MatrixXd levelPoints = V;
		Eigen::MatrixXd levelNormals = normals;
		// Points in coarser level
		Eigen::MatrixXd samplePoints;
		// Neighborhood data structure
		Eigen::MatrixXi neighLevelK = neigh;

		// Compute initial radius
		double densityRatio = std::sqrt(ratio);
		int nLevel1 = int(levelPoints.rows() / ratio);

		// Sampling
		if (samplingStrategy == FPS) {
			samples.push_back(constructFarthestPointSample(V, nLevel1, neigh));
		}
		std::random_device					rd;
		std::default_random_engine			generator(rd());

		// Debug: set number of threads to 1
		//const int NUM_OF_THREADS = omp_get_num_procs();
		const int NUM_OF_THREADS = 1;
		omp_set_num_threads(NUM_OF_THREADS);
		int tid, ntids, ipts, istart, iproc;

		hierarchyTiming["PDS"] = 0.0;
		hierarchyTiming["sampling"] = 0.0;
		hierarchyTiming["cluster"] = 0.0;
		hierarchyTiming["next_neighborhood"] = 0.0;
		hierarchyTiming["next_positions"] = 0.0;
		hierarchyTiming["triangle_finding"] = 0.0;
		hierarchyTiming["triangle_selection"] = 0.0;
		hierarchyTiming["levels"] = DoF.size() - 1.0;
		// For each level
		int k = 0;
		DoF.clear();
		DoF.shrink_to_fit();
		DoF.push_back(levelPoints.rows());
		while (levelPoints.rows() > lowBound && k < 10) {
			double radius = std::cbrt(ratio) * computeAverageEdgeLength(levelPoints, neighLevelK);
			// -- Setting up variables
			// Data structure for neighbors inside level k
			std::vector<std::set<int>> neighborsList;

			// List of triplets to build prolongation operator U
			std::vector<Eigen::Triplet<double>> AllTriplet, UNeighAllTriplet;

			// The nearest coarse point for each fine point and the distances computed
			Eigen::VectorXd D(levelPoints.rows());
			D.setConstant(std::numeric_limits<double>::max());
			nearestSource.push_back(std::vector<size_t>(levelPoints.rows()));

			// -- Sample a subset for level k + 1
			if (verbose) printf("__Constructing Prolongation Operator for level = %d using closest triangle. \n", k);

			// Sample points that will be part of the coarser level with Poisson disk sampling
			if (verbose) std::cout << "Obtaining subset from the finer level\n";
			
			plf::nanotimer samplingTimer;
			samplingTimer.start();

			switch (samplingStrategy) {
				case FASTDISK:
					samples.push_back(fastDiskSample(levelPoints, neighLevelK, radius, D, nearestSource[k]));
					DoF.push_back(samples[k].size());
					break;
				case POISSONDISK:
					DoF.push_back(DoF[k] / ratio);
					samples.push_back(std::vector<int>(0));
					constructPoissonDiskSample(levelPoints, DoF[k + 1], radius, samples[k]);
					break;
				case FPS:
					DoF.push_back(DoF[k] / ratio);
					if (k > 0) {
						samples.push_back(std::vector<int>(DoF[k + 1]));
						std::iota(samples[k].begin(), samples[k].end(), 0);
					}
					break;
				case RANDOM:
					DoF.push_back(DoF[k] / ratio);
					samples.push_back(std::vector<int>(DoF[k]));
					std::iota(samples[k].begin(), samples[k].end(), 0);
					std::shuffle(samples[k].begin(), samples[k].end(), generator);
					samples[k].resize(DoF[k + 1]);
					break;
				case MIS:
					samples.push_back(maximumDeltaIndependentSet(levelPoints, neighLevelK, radius, D, nearestSource[k]));
					DoF.push_back(samples[k].size());
					break;
			}

			if (samples[k].size() < lowBound) {
				nearestSource.pop_back();
				break;
			}

			if (verbose) cout << "Actual number found: " << samples[k].size() << endl;
			DoF[k + 1] = samples[k].size();

			hierarchyTiming["sampling"] += samplingTimer.get_elapsed_ms();

			plf::nanotimer clusterTimer;
			clusterTimer.start();			

			// -- Compute distance from fine points to coarse points and get closest coarse point
			// using distances from MIS if computed before
			constructDijkstraWithCluster(levelPoints, samples[k], neighLevelK, k, D, nearestSource[k]); // Stores result in nearestSource[k]
			hierarchyTiming["cluster"] += clusterTimer.get_elapsed_ms();

			plf::nanotimer nextNeighTimer;
			nextNeighTimer.start();

			// Create neighborhood for the next level
			neighborsList.resize(DoF[k + 1]);
			for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
				for (int j = 0; j < neighLevelK.cols(); ++j) {
					int neighIdx = neighLevelK(fineIdx, j);
					if (neighIdx < 0) break;
					if (nearestSource[k][fineIdx] != nearestSource[k][neighIdx]) {
						neighborsList[nearestSource[k][fineIdx]].insert(nearestSource[k][neighIdx]);
					}
				}
			}

			// Store in homogeneous datastructure
			int maxNeighNum = 0;
			for (int i = 0; i < neighborsList.size(); ++i) {
				if (neighborsList[i].size() > maxNeighNum) {
					maxNeighNum = neighborsList[i].size();
				}
			}
			neighLevelK.resize(DoF[k + 1], maxNeighNum);
			neighLevelK.setConstant(-1);
			for (int i = 0; i < neighborsList.size(); ++i) {
				neighLevelK(i, 0) = i;
				int iCounter = 1;
				for (int node : neighborsList[i]) {
					if (node == i)	continue;
					if (iCounter >= maxNeighNum) break;
					neighLevelK(i, iCounter) = node;
					iCounter++;
				}
			}
			hierarchyTiming["next_neighborhood"] += nextNeighTimer.get_elapsed_ms();

			plf::nanotimer nextPosTimer;
			nextPosTimer.start();
			if (verbose) std::cout << "Setting up the point locations for the next level\n";

			// Setting up the DoF for the next level
			// tempPoints are the centers of the voronoi cells, each row for each voronoi cells
			Eigen::MatrixXd tempPoints(DoF[k + 1], levelPoints.cols());
			tempPoints.setZero();
			if (nested) {
				for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
					tempPoints.row(coarseIdx) = levelPoints.row(samples[k][coarseIdx]);
				}
			} else {
				std::vector<int> clusterSizes(DoF[k + 1]);
				for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
					int coarseIdx = nearestSource[k][fineIdx];
					tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) + levelPoints.row(fineIdx);
					++clusterSizes[coarseIdx];
				}
				for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
					if (clusterSizes[coarseIdx] == 1) {
						tempPoints.row(coarseIdx) = levelPoints.row(samples[k][coarseIdx]);
						for (int neighIdx : neighborsList[coarseIdx]) {
							tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) + levelPoints.row(samples[k][neighIdx]);
						}
						tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) / (neighborsList[coarseIdx].size() + 1.);
					} else {
						tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) / clusterSizes[coarseIdx];
					}
				}
			}
			if (debug) levelV.push_back(tempPoints);
			hierarchyTiming["next_positions"] += nextPosTimer.get_elapsed_ms();

			plf::nanotimer triangleFindingTimer;
			triangleFindingTimer.start();

			// Create triangles for this level based on Voronoi cells
			std::vector<std::vector<int>> tris;
			tris.reserve(DoF[k + 1] * maxNeighNum);
			std::vector<std::vector<int>> connectedTris(DoF[k + 1]);
			std::vector<Eigen::RowVector3d> triNormals;
			triNormals.reserve(DoF[k + 1] * maxNeighNum);
			int triIdx = 0;
			for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
				// Iterate over delaunay triangles
				int v2Idx, v3Idx;
				for (auto it = neighborsList[coarseIdx].begin(); it != neighborsList[coarseIdx].end(); it++) {
					v2Idx = *it;
					// We iterate over the coarse indices in order,
					// so if the neighboring idx is lower then the current coarseIdx,
					// it must have been considered before and be part of a triangle.
					if (v2Idx < coarseIdx) continue;
					for (auto it2 = std::next(it); it2 != neighborsList[coarseIdx].end(); it2++) {
						v3Idx = *it2;
						if (v3Idx < coarseIdx) continue;
						if ((checkVoronoi && neighborsList[v2Idx].find(v3Idx) != neighborsList[v2Idx].end()) || !checkVoronoi) {
							tris.push_back({coarseIdx, v2Idx, v3Idx});
							Eigen::RowVector3d e12 = tempPoints.row(v2Idx) - tempPoints.row(coarseIdx);
							Eigen::RowVector3d e13 = tempPoints.row(v3Idx) - tempPoints.row(coarseIdx);
							triNormals.push_back(e12.cross(e13).normalized());
							connectedTris[coarseIdx].push_back(triIdx);
							connectedTris[v2Idx].push_back(triIdx);
							connectedTris[v3Idx].push_back(triIdx);
							++triIdx;
						}
					}
				}
			}
			tris.shrink_to_fit();
			triNormals.shrink_to_fit();
			if (debug) allTriangles.push_back(tris);	
			hierarchyTiming["triangle_finding"] += triangleFindingTimer.get_elapsed_ms();

			plf::nanotimer triangleSelectionTimer;
			triangleSelectionTimer.start();

			// Create local triangulation on each cluster (centralized at sample i)
			int notrisfound = 0;
			int edgesfound = 0;
			int fallbackCount = 0;
			if (debug) noTriFoundMap.push_back(std::vector<int>(DoF[k]));

			for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
				Eigen::RowVector3d finePoint = levelPoints.row(fineIdx);
				int coarseIdx = nearestSource[k][fineIdx];
				Eigen::RowVector3d coarsePoint = tempPoints.row(coarseIdx);
				std::vector<double> weights;

				if (nested && samples[k][nearestSource[k][fineIdx]] == fineIdx) {
					AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, 1.));
					continue;
				}

				if (neighborsList[coarseIdx].empty()) {
					// If the coarse point has no neighbors,
					// set the weight to 1 for the coarse point.
					// Note: this should not happen.
					AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, 1.));
				} else if (neighborsList[coarseIdx].size() == 1) {
					// If the coarse point only has one neighbor, no triangle can be created.
					// Thus, the weights are distributed w.r.t. the distance to each coarse point.
					int neighIdx = *neighborsList[coarseIdx].begin();
					Eigen::RowVector3d neighPoint = tempPoints.row(neighIdx);

					// get the distance to the two neighboring centroids
					Eigen::RowVector3d e12 = neighPoint - coarsePoint;
					double e12Length = max(e12.norm(), 1e-8);
					double w2 = (levelPoints.row(fineIdx) - coarsePoint).dot(e12.normalized()) / e12Length;
					w2 = min(max(w2, 0.), 1.);
					double w1 = 1 - w2;

					switch (weightingScheme) {
						case BARYCENTRIC:
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, w1));
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, neighIdx, w2));
							break;
						case UNIFORM:
							weights = uniformWeights(2);
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, neighIdx, weights[1]));
							break;
						case INVDIST:
							std::vector<int> endPoints = {coarseIdx, neighIdx};
							weights = inverseDistanceWeights(tempPoints, finePoint, endPoints);
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, neighIdx, weights[1]));
							break;
					}
				} else {
					// Only keep triangle with minimum distance
					double minDistToTriangle = std::numeric_limits<double>::max();
					Eigen::RowVector3d minBary = {0., 0., 0.};
					std::vector<int> minTri;
					bool triFound = false;

					// Values are positive if inside and negative if not
					// Float value represents distance
					std::map<int, float> insideEdge;

					// Iterate over all triangles
					for (int triIdx : connectedTris[coarseIdx]) {
						std::vector<int> tri = tris[triIdx];
						// Make sure that the coarseIdx is in position 0, while keeping orientation
						while(tri[0] != coarseIdx) std::rotate(tri.begin(), tri.begin() + 1, tri.end());

						Eigen::RowVector3d bary = {0., 0., 0.};
						// If the triangle contains the point, the distance is positive, else it's negative
						double distToTriangle = inTriangle(finePoint, tri, triNormals[triIdx], tempPoints, bary, insideEdge);	
						if (distToTriangle >= 0. && distToTriangle < minDistToTriangle) {
							triFound = true;
							minDistToTriangle = distToTriangle;
							minTri = tri;
							minBary = bary;
							break;
						}		
					}

					if (triFound) {
						switch (weightingScheme) {
							case BARYCENTRIC:
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[0], minBary(0)));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[1], minBary(1)));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[2], minBary(2)));
								break;
							case UNIFORM:
								weights = uniformWeights(3);
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[0], weights[0]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[1], weights[1]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[2], weights[2]));
								break;
							case INVDIST:
								weights = inverseDistanceWeights(tempPoints, finePoint, minTri);
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[0], weights[0]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[1], weights[1]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[2], weights[2]));
								break;
						}
					} else {
						bool edgeFound = false;
						double minEdge = std::numeric_limits<double>::max();
						int minEdgeIdx = 0;
						for (const auto& element : insideEdge) {
							const auto& key = element.first;
							const auto& value = element.second;
							if (value >= 0. && value < minEdge) {
								edgeFound = true;
								minEdge = value;
								minEdgeIdx = key;
								break;
							}
						}
						if (edgeFound) {
							++edgesfound;
							Eigen::RowVector3d p2 = tempPoints.row(minEdgeIdx);
							Eigen::RowVector3d e12 = p2 - coarsePoint;
							double e12Length = max(e12.norm(), 1e-8);
							double w2 = (finePoint - coarsePoint).dot(e12.normalized()) / e12Length;
							w2 = min(max(w2, 0.), 1.);
							double w1 = 1. - w2;

							switch (weightingScheme) {
								case BARYCENTRIC:
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, w1));
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minEdgeIdx, w2));
									break;
								case UNIFORM:
									weights = uniformWeights(2);
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minEdgeIdx, weights[1]));
									break;
								case INVDIST:
									std::vector<int> endPointsEdge = {coarseIdx, minEdgeIdx};
									weights = inverseDistanceWeights(tempPoints, finePoint, endPointsEdge);
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minEdgeIdx, weights[1]));
									break;
							}
						} else {
							// Use closest three
							std::vector<int> prolongFrom(3);
							prolongFrom[0] = coarseIdx;

							std::vector<VertexPair> pointsDistances;
							for (int j = 0; j < neighLevelK.cols(); ++j) {
								int neighIdx = neighLevelK(coarseIdx, j);
								if (neighIdx < 0 || neighIdx == coarseIdx) continue;
								VertexPair vp = {neighIdx, (finePoint - tempPoints.row(neighIdx)).norm()};
								pointsDistances.push_back(vp);
							}
							std::sort(pointsDistances.begin(), pointsDistances.end(), std::less<VertexPair>());
							for (int j = 1; j < 3; ++j) {
								prolongFrom[j] = pointsDistances[j - 1].vId;
							}
							std::vector<double> weights = inverseDistanceWeights(tempPoints, finePoint, prolongFrom);
							for (int j = 0; j < prolongFrom.size(); j++) {
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, prolongFrom[j], weights[j]));
							}
							++fallbackCount;
						}
						++notrisfound;
					}
				}
			}
			if (verbose) cout << "Number of triangles unfound: " << notrisfound << endl;
			if (verbose) cout << "Number of edges found: " << edgesfound << endl;
			if (verbose) cout << "Percentage of fallback: " << (double) fallbackCount / (double) DoF[k] * 100 << endl;
			hierarchyTiming["triangle_selection"] += triangleSelectionTimer.get_elapsed_ms();

			levelPoints = tempPoints;

			Eigen::SparseMatrix<double> ULevel;
			ULevel.resize(DoF[k], DoF[k + 1]);
			ULevel.setFromTriplets(AllTriplet.begin(), AllTriplet.end());
			U.push_back(ULevel);
			AllTriplet.clear();
			AllTriplet.shrink_to_fit();
			++k;
		}
	}

	double MultigridSolver::inTriangle(const Eigen::RowVector3d& p, const std::vector<int>& tri, const Eigen::RowVector3d& triNormal, const Eigen::MatrixXd& pos, Eigen::RowVector3d& bary, std::map<int, float>& insideEdge)
	{
		Eigen::RowVector3d v1, v2, v3;
		v1 = pos.row(tri[0]);
		v2 = pos.row(tri[1]);
		v3 = pos.row(tri[2]);
		Eigen::RowVector3d v1ToP = p - v1;
		Eigen::RowVector3d e12 = v2 - v1;
		Eigen::RowVector3d e13 = v3 - v1;

		double distToTriangle = (p - v1).dot(triNormal);
		Eigen::RowVector3d pProjected = p - distToTriangle * triNormal;

		double doubleArea = (v2 - v1).cross(v3 - v1).dot(triNormal);
		bary(0) = (v3 - v2).cross(pProjected - v2).dot(triNormal) / doubleArea;
		bary(1) = (v1 - v3).cross(pProjected - v3).dot(triNormal) / doubleArea;
		bary(2) = 1. - bary(0) - bary(1);
		
		if (insideEdge.find(tri[1]) == insideEdge.end()) {
			insideEdge[tri[1]] = ((v1ToP) - ((v1ToP).dot(e12) * (e12))).norm();
		}
		if (insideEdge.find(tri[2]) == insideEdge.end()) {
			insideEdge[tri[2]] = ((v1ToP) - ((v1ToP).dot(e13) * (e13))).norm();
		}
		if (bary(0) < 0. || bary(1) < 0.) {
			insideEdge[tri[1]] = -1.;
		}
		if (bary(0) < 0. || bary(2) < 0.) {
			insideEdge[tri[2]] = -1.;
		}

		if (bary(0) >= 0. && bary(1) >= 0. && bary(2) >= 0.) {
			return abs(distToTriangle);
		}

		return -1.;
	}

	std::vector<double> MultigridSolver::uniformWeights(const int& n_points) {
		std::vector<double> weights(n_points);
		std::fill(weights.begin(), weights.end(), 1. / n_points);
		return weights;
	}

	std::vector<double> MultigridSolver::inverseDistanceWeights(const Eigen::MatrixXd& pos, const Eigen::RowVector3d& p, const std::vector<int>& edges) {
		double sumWeight = 0.;
		std::vector<double> weights(edges.size());
		for (int j = 0; j < edges.size(); ++j) {
			weights[j] = 1. / max(1e-8, (p - pos.row(edges[j])).norm());
			sumWeight += weights[j];
		}
		for (int j = 0; j < weights.size(); ++j) {
			weights[j] = weights[j] / sumWeight;
		}
		return weights;
	}

	void MultigridSolver::constructProlongationSIG06()
	{
		// Create prolongation operator for level k+1 to k
		// Points in current level
		Eigen::MatrixXd levelPoints = V;
		// Neighborhood data structure
		Eigen::MatrixXi neighLevelK = neigh;
		// Initial radius graph

		// Sampling
		std::random_device					rd;
		std::default_random_engine			generator(rd());

		// Debug: set number of threads to 1
		//const int NUM_OF_THREADS = omp_get_num_procs();
		const int NUM_OF_THREADS = 1;
		omp_set_num_threads(NUM_OF_THREADS);
		int tid, ntids, ipts, istart, iproc;

		hierarchyTiming["PDS"] = 0.0;
		hierarchyTiming["sampling"] = 0.0;
		hierarchyTiming["cluster"] = 0.0;
		hierarchyTiming["next_neighborhood"] = 0.0;
		hierarchyTiming["next_positions"] = 0.0;
		hierarchyTiming["triangulation"] = 0.0;
		hierarchyTiming["levels"] = DoF.size() - 1.0;
		// For each level
		int k = 0;
		DoF.clear();
		DoF.shrink_to_fit();
		DoF.push_back(levelPoints.rows());
		while (levelPoints.rows() > lowBound && k < 10) {
			double radius = std::cbrt(5.) * computeAverageEdgeLength(levelPoints, neighLevelK);
			// -- Setting up variables
			// Data structure for neighbors inside level k
			std::vector<std::set<int>> neighborsList;

			// List of triplets to build prolongation operator U
			std::vector<Eigen::Triplet<double>> AllTriplet, UNeighAllTriplet;

			// -- Sample a subset for level k + 1
			if (verbose) printf("__Constructing Prolongation Operator for level = %d using closest triangle. \n", k);

			// Sample points that will be part of the coarser level with Poisson disk sampling
			if (verbose) std::cout << "Obtaining subset from the finer level\n";
			
			plf::nanotimer pdsTimer;
			pdsTimer.start();
			
			// Create a maximum del-independent set with target fraction 0.2,
			// as described in SIG06 paper
			samples.push_back(maximumDeltaIndependentSet(levelPoints, neighLevelK, radius));
			if (samples[k].size() < lowBound || (k > 1 && (double) samples[k].size() / (double) samples[k - 1].size() > 0.9)) break;

			if (verbose) cout << "Actual number found: " << samples[k].size() << endl;

			// In case the number of samples is not exactly the same, update:
			DoF.push_back(samples[k].size());
			neighborsList.resize(DoF[k + 1]);
			// Store map from original point indices to sample indices
			PointsToSampleMap.resize(levelPoints.rows(), -1);
			std::fill(PointsToSampleMap.begin(), PointsToSampleMap.end(), -1);
			for (int i = 0; i < samples[k].size(); ++i) {
				PointsToSampleMap[samples[k][i]] = i;
			}
			PointsToSampleMaps.push_back(PointsToSampleMap);

			hierarchyTiming["PDS"] += pdsTimer.get_elapsed_ms();

			plf::nanotimer nextNeighTimer;
			nextNeighTimer.start();

			// Create coarser graph by going through 2-ring
			// And store new coarse node positions in tempPoints
			Eigen::MatrixXd tempPoints(DoF[k + 1], levelPoints.cols());
			for (int i = 0; i < DoF[k + 1]; ++i) {
				tempPoints.row(i) = levelPoints.row(samples[k][i]);
				for (int j = 0; j < neighLevelK.cols(); ++j) {
					int neigh1Ring = neighLevelK(samples[k][i], j);
					if (neigh1Ring < 0) break;
					if (PointsToSampleMap[neigh1Ring] < DoF[k + 1] && PointsToSampleMap[neigh1Ring] != i) neighborsList[i].insert(PointsToSampleMap[neigh1Ring]);
					for (int j2 = 0; j2 < neighLevelK.cols(); ++j2) {
						int neigh2Ring = neighLevelK(neigh1Ring, j2);
						if (neigh2Ring < 0) break;
						if (PointsToSampleMap[neigh2Ring] < DoF[k + 1] && PointsToSampleMap[neigh2Ring] != i) neighborsList[i].insert(PointsToSampleMap[neigh2Ring]);
					}
				}
			}
			levelV.push_back(tempPoints);

			hierarchyTiming["next_neighborhood"] += nextNeighTimer.get_elapsed_ms();

			plf::nanotimer triangulationTimer;
			triangulationTimer.start();

			// Assign weights to neighbors in original graph
			for (int i = 0; i < DoF[k]; ++i) {
				if (PointsToSampleMap[i] < DoF[k + 1]) {
					AllTriplet.push_back(Eigen::Triplet<double>(i, PointsToSampleMap[i], 1.));
					continue;
				}
				std::vector<int> coarseNeighs;
				coarseNeighs.reserve(neighLevelK.cols());
				std::vector<double> coarseWeights;
				coarseWeights.reserve(neighLevelK.cols());
				double sumWeight = 0;
				for (int j = 0; j < neighLevelK.cols(); ++j) {
					int neighIdx = neighLevelK(i, j);
					if (neighIdx < 0) break;
					if (PointsToSampleMap[neighIdx] < DoF[k + 1]) {
						coarseNeighs.push_back(PointsToSampleMap[neighIdx]);
						double weight = 1. / max(1e-8, (levelPoints.row(i) - levelPoints.row(neighIdx)).norm());
						sumWeight += weight;
						coarseWeights.push_back(weight);
					}
				}
				coarseNeighs.shrink_to_fit();
				coarseWeights.shrink_to_fit();
				double numCoarse = (double) coarseNeighs.size();
				for (int j = 0; j < coarseNeighs.size(); ++j) {
					AllTriplet.push_back(Eigen::Triplet<double>(i, coarseNeighs[j], coarseWeights[j] / sumWeight));
				}
			}
			hierarchyTiming["triangulation"] += triangulationTimer.get_elapsed_ms();

			// Create neighborhood for the next level
			int maxNeighNum = 0;
			int totalE = 0;
			for (int i = 0; i < neighborsList.size(); ++i) {
				totalE += neighborsList[i].size();
				if (neighborsList[i].size() > maxNeighNum) {
					maxNeighNum = neighborsList[i].size();
				}
			}

			Eigen::MatrixXi levelEdges;
			levelEdges.resize(totalE, 2);
			levelEdges.setConstant(0);
			int countEdges = 0;
			neighLevelK.resize(DoF[k + 1], maxNeighNum);
			neighLevelK.setConstant(-1);
			for (int i = 0; i < neighborsList.size(); ++i) {
				int iCounter = 0;
				for (int node : neighborsList[i]) {
					levelEdges(countEdges, 0) = i;
					levelEdges(countEdges, 1) = node;
					++countEdges;
					if (node == i)	continue;
					if (iCounter >= maxNeighNum) break;
					neighLevelK(i, iCounter) = node;
					iCounter++;
				}
			}
			levelE.push_back(levelEdges);

			levelPoints = tempPoints;

			Eigen::SparseMatrix<double> ULevel;
			ULevel.resize(DoF[k], DoF[k + 1]);
			ULevel.setFromTriplets(AllTriplet.begin(), AllTriplet.end());
			U.push_back(ULevel);
			AllTriplet.clear();
			AllTriplet.shrink_to_fit();
			++k;
		}
	}

	double MultigridSolver::computeAverageEdgeLength(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& neigh) {
		double sumLength = 0;
		int nEdges = 0;
		for (int i = 0; i < pos.rows(); ++i) {
			Eigen::Vector3d p1 = pos.row(i);
			for (int j = 0; j < neigh.cols(); ++j) {
				if (neigh(i, j) < 0) continue;
				Eigen::Vector3d p2 = pos.row(neigh(i, j));
				double dist = (p1 - p2).norm();
				if (dist > 0) {
					sumLength += dist;
					++nEdges;
				}
			}
		}
		return sumLength / (double) nEdges;
	}

	//!< Construct a regular Poisson Disk Sampling 
	void MultigridSolver::constructPoissonDiskSample(const Eigen::MatrixXd& points, const int sampleSize, const float radius, std::vector<int>& PDSSamples)
	{
		int nBox;
		if (sampleSize < 1 || sampleSize > points.rows()) {
			std::cout << "Error! Sample cannot be 0 or larger than the number of vertices." << std::endl;
			return;
		}
		else if (sampleSize <= 300) {
			nBox = 8;
		}
		else if (sampleSize > 300 && sampleSize <= 500) {
			nBox = 10;
		}
		else if (sampleSize > 500 && sampleSize <= 2000) {
			nBox = 13;
		}
		else if (sampleSize > 2000 && sampleSize <= 3000) {
			nBox = 16;
		}
		else if (sampleSize > 3000 && sampleSize <= 6000) {
			nBox = 20;
		}
		else if (sampleSize > 6000 && sampleSize <= 10000) {
			nBox = 30;
		}
		else if (sampleSize > 10000 && sampleSize <= 30000) {
			nBox = 40;
		}
		else if (sampleSize > 30000 && sampleSize <= 75000) {
			nBox = 75;
		}
		else if (sampleSize > 75000 && sampleSize <= 100000) {
			nBox = 100;
		}
		else {
			nBox = 150;
		}

		/* Initialization of computing boundary of each box	 */
		Eigen::Vector3d minV(points.col(0).minCoeff(), points.col(1).minCoeff(), points.col(2).minCoeff());
		Eigen::Vector3d maxV(points.col(0).maxCoeff(), points.col(1).maxCoeff(), points.col(2).maxCoeff());
		Eigen::Vector3d range, boxDist;
		Eigen::Vector3i nn;

		range = maxV - minV;
		double minRange = range.minCoeff();

		for (int i = 0; i < 3; i++) {
			boxDist(i) = minRange / (double)nBox;
			nn(i) = (int)ceil(range(i) / boxDist(i));
		}

		set<int>		locSample;
		const int		SamplePerBox = 10;

		/* Container for boxes */
		priority_queue<BoxStruct, std::vector<BoxStruct>, std::greater<BoxStruct>>	UnvisitedBoxQueue;
		list<int>																	VisitedEmptyBoxList;

		/* Constructing tensor (n x n x n x (box-members)) for sample candidates */
		list<int>							v1;
		vector<list<int>>					v2;
		vector<vector<list<int>>>			v3;
		vector<vector<vector<list<int>>>>	v4;

		/* Instatiating the tensor with empty elements */
		v2.reserve(nn(2));	v3.reserve(nn(1));	v4.reserve(nn(0));
		for (int i = 0; i < nn(2); i++) { v2.push_back(v1); }
		for (int i = 0; i < nn(1); i++) { v3.push_back(v2); }
		for (int i = 0; i < nn(0); i++) { v4.push_back(v3); }

		/* Put each vertex to its corresponding box */
		for (int vId = 0; vId < points.rows(); vId++) {
			int xId, yId, zId;
			xId = (int)abs(floor((points(vId, 0) - minV(0)) / boxDist(0)));
			yId = (int)abs(floor((points(vId, 1) - minV(1)) / boxDist(1)));
			zId = (int)abs(floor((points(vId, 2) - minV(2)) / boxDist(2)));
			if (xId > (nn(0) - 1)) xId = nn(0) - 1;
			if (yId > (nn(1) - 1)) yId = nn(1) - 1;
			if (zId > (nn(2) - 1)) zId = nn(2) - 1;

			v4[xId][yId][zId].push_back(vId);
		}

		// Randomizer for large number (LONG)
		std::random_device					rd;
		std::default_random_engine			generator(rd());
		int									trimCounter = 0;

		/* FOR TESTING ONLY*/
		int poolCounter = 0;
		/* END OF TESTING */

		// Populating Boxes with 10 random vertices
		for (int iIdx = 0; iIdx < nn(0); iIdx++) {
			for (int jIdx = 0; jIdx < nn(1); jIdx++) {
				for (int kIdx = 0; kIdx < nn(2); kIdx++) {
					Eigen::Vector3d		boxMin(minV(0) + iIdx * boxDist(0), minV(1) + jIdx * boxDist(1), minV(2) + kIdx * boxDist(2));

					/* IF there are more samples on a box then required, select randomly 10 of them; */
					int boxSize = v4[iIdx][jIdx][kIdx].size();
					if (boxSize == 0) {
						// ID(a,b,c) = a*(y*z) + b*(z) + c
						VisitedEmptyBoxList.push_back(iIdx * nn(2) * nn(1) + jIdx * nn(2) + kIdx);
					}
					else
					{	// Work on NON-EMPTY boxes only
						double randDist = (double)(rand() % 100) / (double)10;
						BoxStruct curBox{ iIdx * nn(2) * nn(1) + jIdx * nn(2) + kIdx, randDist };
						UnvisitedBoxQueue.push(curBox);

						/* If there are more elements than required, randomly discard the rest */
						if (boxSize > SamplePerBox) {
							set<int> vertexToDel;
							srand(time(NULL));
							do {
								std::uniform_int_distribution<long> distribution(0, boxSize);
								int vSel = distribution(generator);
								vertexToDel.insert(vertexToDel.end(), vSel);
							} while (vertexToDel.size() <= (boxSize - SamplePerBox));

							int l = 0;
							for (list<int>::iterator it = v4[iIdx][jIdx][kIdx].begin(); it != v4[iIdx][jIdx][kIdx].end(); ) {
								if (vertexToDel.find(l) != vertexToDel.end()) {
									it = v4[iIdx][jIdx][kIdx].erase(it);
								}
								else {
									++it;
								}
								l++;
							}
							trimCounter++;
						} // End of vertex(-ices) deletion
					}
					/* FOR TESTING ONLY*/
					poolCounter += v4[iIdx][jIdx][kIdx].size();
					/* END OF TESTING */
				}
			}
		} // End of Populating the boxes with samples POOL

		/* Getting samples from POOL */
		while (VisitedEmptyBoxList.size() < nn(0) * nn(1) * nn(2)) {
			/* Randomly select a certain box */
			BoxStruct box1 = UnvisitedBoxQueue.top();
			UnvisitedBoxQueue.pop();

			/* Determining in which box is this point located */
			int b_xId, b_yId, b_zId;
			b_xId = round(box1.id / (nn(1) * nn(2)));
			b_yId = round((box1.id - b_xId * nn(1) * nn(2)) / nn(2));
			b_zId = box1.id % nn(2);
			if (v4[b_xId][b_yId][b_zId].size() < 1) continue;

			/* Picking a random sample from a box */
			srand(time(NULL));
			int randV1 = rand() % v4[b_xId][b_yId][b_zId].size();
			int boxSample;
			int locCounter = 0;
			for (list<int>::iterator it = v4[b_xId][b_yId][b_zId].begin(); it != v4[b_xId][b_yId][b_zId].end(); ++it) {
				if (randV1 == locCounter) boxSample = *it;
				locCounter++;
			}
			locSample.insert(locSample.end(), boxSample);

			/* Iterating through neigboring boxes */
			for (int i = b_xId - 1; i <= b_xId + 1; i++) {
				if (i < 0 || i >(nn(0) - 1)) continue;
				for (int j = b_yId - 1; j <= b_yId + 1; j++) {
					if (j < 0 || j >(nn(1) - 1)) continue;
					for (int k = b_zId - 1; k <= b_zId + 1; k++) {
						if (k < 0 || k >(nn(2) - 1)) continue;
						if (v4[i][j][k].size() > 0) {
							/* Removing samples closer than a certain radius */
							for (list<int>::iterator it = v4[i][j][k].begin(); it != v4[i][j][k].end(); ) {
								double distAB = (points.row(boxSample) - points.row(*it)).norm();
								//VtoVDist(V.row(boxSample), V.row(*it), distAB);

								if (distAB <= radius) {
									it = v4[i][j][k].erase(it);
								}
								else {
									++it;
								}
							} /* End of Removing samples closer than a certain radius */

								/* Put the empty boxes to its list */
							if (v4[i][j][k].size() < 1) {
								VisitedEmptyBoxList.push_back(i * nn(1) * nn(2) + j * nn(2) + k);
							}
							else
							{
								/* If the original box is not empty, put it back to queue with new priority */
								if (i == b_xId && j == b_yId && k == b_zId) {
									double nRandDist = (double)(rand() % 100) / (double)10;
									// ID = a*(y*z) + b*(z) + c.
									BoxStruct curBox{ i * nn(1) * nn(2) + j * nn(2) + k, nRandDist };
									UnvisitedBoxQueue.push(curBox);
								}
							}
						}
					} /* End of iteration on each box */
				}
			} /* End of Iterating through neigboring boxes */
		}
		/* END of Iterating to get All samples */

		/* Passing sample to caller */
		PDSSamples.resize(locSample.size());

		int sId = 0;
		for (int samp : locSample) {
			PDSSamples[sId++] = samp;
		}
	}

	std::vector<int> MultigridSolver::maximumDeltaIndependentSet(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& edges, const double& radius) {
		std::vector<bool> visited(edges.rows());
		std::vector<int> selection;
		for (int i = 0; i < edges.rows(); ++i) {
			if (!visited[i]) {
				selection.push_back(i);
				for (int j = 0; j < edges.cols(); ++j) {
					int neighIdx = edges(i, j);
					if (neighIdx < 0) break;
					double dist = (pos.row(i) - pos.row(neighIdx)).norm();
					if (dist < radius) {
						visited[neighIdx] = true;
					}
				}
			}
		}
		return selection;
	}

	std::vector<int> MultigridSolver::maximumDeltaIndependentSet(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& edges, const double& radius, Eigen::VectorXd& D, std::vector<size_t>& nearestSourceK) {
		std::vector<bool> visited(edges.rows());
		std::vector<int> selection;
		int sampleIdx = 0;
		for (int i = 0; i < edges.rows(); ++i) {
			if (!visited[i]) {
				selection.push_back(i);
				nearestSourceK[i] = sampleIdx;
				for (int j = 0; j < edges.cols(); ++j) {
					int neighIdx = edges(i, j);
					if (neighIdx < 0) break;
					double dist = (pos.row(i) - pos.row(neighIdx)).norm();
					if (dist < radius) {
						visited[neighIdx] = true;
						if (dist < D(neighIdx)) {
							D(neighIdx) = dist;
							nearestSourceK[neighIdx] = sampleIdx;
						}
					}
				}
				++sampleIdx;
			}
		}
		return selection;
	}

	std::vector<int> MultigridSolver::fastDiskSample(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& edges, const double& radius, Eigen::VectorXd& D, std::vector<size_t>& nearestSourceK) {
		std::vector<bool> visited(edges.rows());
		std::vector<int> selection;
		int sampleIdx = 0;
		for (int i = 0; i < edges.rows(); ++i) {
			if (!visited[i]) {
				selection.push_back(i);
				nearestSourceK[i] = sampleIdx;
				for (int j = 0; j < edges.cols(); ++j) {
					int neighIdx = edges(i, j);
					if (neighIdx < 0) break;
					double dist = (pos.row(i) - pos.row(neighIdx)).norm();
					if (dist < radius) {
						visited[neighIdx] = true;
						if (dist < D(neighIdx)) {
							D(neighIdx) = dist;
							nearestSourceK[neighIdx] = sampleIdx;
						}
						if (!sig06) {
							for (int j2 = 0; j2 < edges.cols(); ++j2) {
								int neighIdx2 = edges(neighIdx, j2);
								if (neighIdx2 < 0) break;
								double dist2 = dist + (pos.row(neighIdx) - pos.row(neighIdx2)).norm();
								if (dist2 < radius) {
									visited[neighIdx2] = true;
									if (dist2 < D(neighIdx2)) {
										D(neighIdx2) = dist2;
										nearestSourceK[neighIdx2] = sampleIdx;
									}
								}
							}
						}
					}
				}
				++sampleIdx;
			}
		}
		return selection;
	}

	void MultigridSolver::constructDijkstraWithCluster(const Eigen::MatrixXd& points, const std::vector<int>& source, const Eigen::MatrixXi& neigh, int k, Eigen::VectorXd& D, std::vector<size_t>& nearestSourceK)
	{
		std::priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistanceQueue;
		if (nearestSourceK.empty())	nearestSourceK.resize(points.rows(), source[0]);

		for (int i = 0; i < source.size(); ++i) {
			D(source[i]) = 0.0;
			VertexPair vp{ source[i], D(source[i]) };
			DistanceQueue.push(vp);
			nearestSourceK[source[i]] = i;
		}

		int curSource;
		while (!DistanceQueue.empty()) {
			VertexPair vp1 = DistanceQueue.top();
			curSource = nearestSourceK[vp1.vId];
			Eigen::RowVector3d vertex1 = points.row(vp1.vId);
			DistanceQueue.pop();

			//for (int vNeigh : neigh.row(vp1.vId)) {
			for (int i = 0; i < neigh.cols(); ++i) {
				int vNeigh = neigh(vp1.vId, i);

				if (vNeigh >= 0) {
					double dist, distTemp;
					Eigen::RowVector3d vertex2 = points.row(vNeigh);
					dist = (vertex2 - vertex1).norm();
					distTemp = vp1.distance + dist;
					if (distTemp < D(vNeigh)) {
						// Assign a new distance
						D(vNeigh) = distTemp;
						VertexPair v2{ vNeigh, distTemp };
						DistanceQueue.push(v2);


						// Assign the nearest source to a certain point
						nearestSourceK[vNeigh] = curSource;
					}
				}
			}
		}
	}


	double  MultigridSolver::multiGridVCycleGS(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& b, Eigen::MatrixXd& x, int k, bool isDebug)
	{
		double tolPre = 1e-15, tolPost = 1e-15;
		//[1] Smoothing
		GaussSeidelSmoother(A, b, x, preIters, tolPre, isDebug);

		//[2] Residual
		Eigen::MatrixXd res = b - A * x;

		//[3] Restriction
		Eigen::MatrixXd resRest = U[k].transpose() * res;

		//[4] Recursion
		Eigen::MatrixXd eps;
		eps.setZero(resRest.rows(), resRest.cols());
		if (k == U.size() - 1) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridVCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[5] Prolongation
		x = x + U[k] * eps;

		//[6] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		return 0.0;
	}

	//!< Perform Multigrid using F-Cycle and Gauss-Seidel smoothing
	double  MultigridSolver::multiGridFCycleGS(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& b, Eigen::MatrixXd& x, int k, bool isDebug)
	{
		double tolPre = 1e-15, tolPost = 1e-15;
		//[1] Smoothing
		GaussSeidelSmoother(A, b, x, preIters, tolPre, isDebug);

		//[2] Residual
		Eigen::MatrixXd res = b - A * x;

		//[3] Restriction
		Eigen::MatrixXd resRest = U[k].transpose() * res;

		//[4] Recursion
		Eigen::MatrixXd eps;
		eps.setZero(resRest.rows(), resRest.cols());
		if (k == U.size() - 1) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridFCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[5] Prolongation
		x = x + U[k] * eps;

		//[6] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		//[7] Residual
		res = b - A * x;

		//[8] Restriction
		resRest = U[k].transpose() * res;

		//[9] Recursion
		if (k == DoF.size() - 2) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridVCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[10] Prolongation
		x = x + U[k] * eps;

		//[11] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		return 0.0;
	}

	//!< Perform Multigrid using W-Cycle and Gauss-Seidel smoothing
	double  MultigridSolver::multiGridWCycleGS(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& b, Eigen::MatrixXd& x, int k, bool isDebug)
	{
		double tolPre = 1e-15, tolPost = 1e-15;
		//[1] Smoothing
		GaussSeidelSmoother(A, b, x, preIters, tolPre, isDebug);

		//[2] Residual
		Eigen::MatrixXd res = b - A * x;

		//[3] Restriction
		Eigen::MatrixXd resRest = U[k].transpose() * res;

		//[4] Recursion
		Eigen::MatrixXd eps;
		eps.setZero(resRest.rows(), resRest.cols());
		if (k == U.size() - 1) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridWCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[5] Prolongation
		x = x + U[k] * eps;

		//[6] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		//[7] Residual
		res = b - A * x;

		//[8] Restriction
		resRest = U[k].transpose() * res;

		//[9] Recursion
		if (k == DoF.size() - 2) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridWCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[10] Prolongation
		x = x + U[k] * eps;

		//[11] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		return 0.0;
	}

	void MultigridSolver::GaussSeidelSmoother(Eigen::SparseMatrix<double>& LHS, Eigen::MatrixXd& rhs,
		Eigen::MatrixXd& x, int maxIter, double tol, bool isDebug)
	{
		int dim = x.cols();
		if (dim == 1) {
			for (int i = 0; i < maxIter; i++) {
				for (int k = 0; k < LHS.outerSize(); ++k) {
					double sum = 0.0;
					for (Eigen::SparseMatrix<double>::InnerIterator it(LHS, k); it; ++it) {
						if (it.row() != k) {
							sum += it.value() * x(it.row());
						}
					}
					x(k) = (rhs(k) - sum) / LHS.coeffRef(k, k);
				}
			}
		}
		else {
			for (int i = 0; i < maxIter; i++) {
				for (int iCol = 0; iCol < dim; ++iCol) {
					for (int k = 0; k < LHS.outerSize(); ++k) {
						double sum = 0.0;
						for (Eigen::SparseMatrix<double>::InnerIterator it(LHS, k); it; ++it) {
							if (it.row() != k) {
								sum += it.value() * x(it.row(), iCol);
							}
						}
						x(k, iCol) = (rhs(k, iCol) - sum) / LHS.coeffRef(k, k);
					}
				}
			}
		}
	}

	double MultigridSolver::residualCheck(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& b, const Eigen::MatrixXd& x, int type)
	{
		// 0: rel. norm (||Ax-b||/||b||)		1: L2 M^-1 (||Ax-b||_{M-1}/||b||_{M-1})		2: L2 M (||Ax-b||{M}/||b||_{M})		3: Abs (Ax-b).norm()
		double residue;
		switch (type)
		{
		case 0:
		{
			Eigen::VectorXd absResidue(b.cols());
			for (int i = 0; i < b.cols(); ++i) {
				absResidue(i) = (A * x.col(i) - b.col(i)).norm() / b.col(i).norm();
			}
			residue = absResidue.maxCoeff();
			break;
		}
		case 1:
		{
			Eigen::VectorXd absResidue(b.cols());
			Eigen::VectorXd residual;
			for (int i = 0; i < b.cols(); ++i) {
				residual = A * x.col(i) - b.col(i);
				double n1 = residual.transpose() * Minv * residual;
				double n2 = b.col(i).transpose() * Minv * b.col(i);
				absResidue(i) = sqrt(n1 / n2);
			}
			residue = absResidue.maxCoeff();
			break;
		}
		case 2:
		{
			Eigen::VectorXd absResidue(b.cols());
			Eigen::VectorXd residual;
			for (int i = 0; i < b.cols(); ++i) {
				residual = A * x.col(i) - b.col(i);
				double n1 = residual.transpose() * M * residual;
				double n2 = b.col(i).transpose() * M * b.col(i);
				absResidue(i) = sqrt(n1 / n2);
			}
			residue = absResidue.maxCoeff();
			break;
		}
		case 3:
			residue = (A * x - b).norm();
			break;
		default:
			break;
		}

		return residue;
	}

	void MultigridSolver::solve(Eigen::SparseMatrix<double>& LHS, Eigen::MatrixXd& rhs, Eigen::MatrixXd& x, int solverType)
	{
		std::chrono::high_resolution_clock::time_point t0, t1, t2;
		std::chrono::duration<double> duration1, duration2;
		convergence.reserve(maxIter);

		Eigen::MatrixXd mx = rhs;

		if (solverType == 0) {
			if (verbose) std::cout << "Solve the linear system using direct solver";

			plf::nanotimer directFactorTimer, directSolveTimer;
			if (factorizationType == 0) {
				if (verbose) std::cout << "(LLT)\n";
				directFactorTimer.start();
				Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> sparseSolver;
				sparseSolver.compute(LHS);
				solverTiming["direct_factor"] = directFactorTimer.get_elapsed_ms();

				directSolveTimer.start();
				x = sparseSolver.solve(mx);
				solverTiming["direct_solve"] = directSolveTimer.get_elapsed_ms();
				solverTiming["direct_residual"] = residualCheck(LHS, mx, x, stoppingCriteria);

				if (verbose) std::cout << "Factor: " << solverTiming["direct_factor"] << std::endl;
				if (verbose) std::cout << "Solve: " << solverTiming["direct_solve"] << std::endl;
			}
			else {
				if (verbose) std::cout << "(LDLT)\n";

				directFactorTimer.start();
				Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver;
				sparseSolver.compute(LHS);
				solverTiming["direct_factor"] = directFactorTimer.get_elapsed_ms();

				directSolveTimer.start();
				x = sparseSolver.solve(mx);
				solverTiming["direct_solve"] = directSolveTimer.get_elapsed_ms();

				if (verbose) std::cout << "Factor: " << solverTiming["direct_factor"] << std::endl;
				if (verbose) std::cout << "Solve: " << solverTiming["direct_solve"] << std::endl;
			}
		}
		else if (solverType == 1) {
			#ifdef PARDISO_ENABLED
			if (verbose) std::cout << "Solve the linear system using Intel MKL's PARDISO direct solver ";

			plf::nanotimer pardisoFactorTimer, pardisoSolveTimer;

			if (factorizationType == 0) {
				if (verbose) std::cout << "(LLT) \n";

				Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver;
				if (verbose) std::cout << "Factorize \n";
				//sparseSolver.compute(LHS);
				pardisoFactorTimer.start();
				sparseSolver.analyzePattern(LHS);
				sparseSolver.factorize(LHS);
				solverTiming["pardiso_factor"] = pardisoFactorTimer.get_elapsed_ms();

				if (verbose) std::cout << "Solve \n";
				pardisoSolveTimer.start();
				x = sparseSolver.solve(mx);
				solverTiming["pardiso_solve"] = pardisoSolveTimer.get_elapsed_ms();

				if (verbose) std::cout << "Factor: " << solverTiming["pardiso_factor"] << std::endl;
				if (verbose) std::cout << "Solve: " << solverTiming["pardiso_solve"] << std::endl;
			}
			else {
				if (verbose) std::cout << "(LDLT) \n";

				Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver;

				pardisoFactorTimer.start();
				sparseSolver.compute(LHS);
				solverTiming["pardiso_factor"] = pardisoFactorTimer.get_elapsed_ms();

				pardisoSolveTimer.start();
				x = sparseSolver.solve(mx);
				solverTiming["pardiso_solve"] = pardisoSolveTimer.get_elapsed_ms();

				if (verbose) std::cout << "Factor: " << solverTiming["pardiso_factor"] << std::endl;
				if (verbose) std::cout << "Solve: " << solverTiming["pardiso_solve"] << std::endl;
			}
			#else // PARDISO_ENSABLED
			if (verbose) std::cout << "Cannot use Pardiso, since the Intel OpenMKL library was not found.";	
			#endif // PARDISO_DISABLED
		}
		else if (solverType == 2) {
			if (verbose) std::cout << "Solve the linear system using OUR Multigrid solver \n";

			double residue = std::numeric_limits<double>::max();
			Eigen::VectorXd rhs0;
			Eigen::VectorXd rhs1;
			Eigen::VectorXd rhs2;
			rhs0 = mx.col(0);
			if (mx.cols() > 1) {
				rhs1 = mx.col(1);
				rhs2 = mx.col(2);
			}

			if (isSmootherGaussSeidel)
			{
				plf::nanotimer solverTotalTimer, reductionTimer;
				solverTotalTimer.start(); reductionTimer.start();

				if (verbose) cout << "Reducing system" << endl;

				Abar.resize(U.size() + 1);
				//Mbar.resize(DoF.size());
				Abar[1] = U[0].transpose() * LHS * U[0];
				for (int k = 2; k < U.size() + 1; ++k) {
					Abar[k] = U[k - 1].transpose() * Abar[k - 1] * U[k - 1];
				}

				solverTiming["reduction"] = reductionTimer.get_elapsed_ms();

				plf::nanotimer coarseSolveTimer;
				coarseSolveTimer.start();

				if (verbose) cout << "Solving reduced system" << endl;

				coarsestSolver.compute(Abar[U.size()]);

				solverTiming["coarsest_solve"] = coarseSolveTimer.get_elapsed_ms();

				plf::nanotimer cyclesTimer;
				cyclesTimer.start();
				int iterCount = 0;
				if (cycleType == 0) {
					if (verbose) std::cout << "V-CYCLE \n";
					Eigen::Vector3d relResidue;
					do {
						residue = multiGridVCycleGS(LHS, mx, x, 0, false);
						residue = residualCheck(LHS, mx, x, stoppingCriteria);
						convergence.push_back({cyclesTimer.get_elapsed_ms(), residue});
						if (verbose) printf("%d,%f,%.14f \n", ++iterCount, cyclesTimer.get_elapsed_ms(), residue);
						else ++iterCount;
					} while ((residue > accuracy) && (iterCount < maxIter));
					convergence.shrink_to_fit();
				}
				else if (cycleType == 1) {
					if (verbose) std::cout << "F-CYCLE \n";
					Eigen::Vector3d relResidue;
					do {
						residue = multiGridFCycleGS(LHS, mx, x, 0, false);
						residue = residualCheck(LHS, mx, x, stoppingCriteria);
						if (verbose) printf("Iteration %d: residue=%.10f \n", ++iterCount, residue);
						else ++iterCount;
					} while ((residue > accuracy) && (iterCount < maxIter));
				}
				else if (cycleType == 2) {
					if (verbose) std::cout << "W-CYCLE \n";
					Eigen::Vector3d relResidue;
					do {
						residue = multiGridWCycleGS(LHS, mx, x, 0, false);
						residue = residualCheck(LHS, mx, x, stoppingCriteria);
						if (verbose) printf("Iteration %d: residue=%.10f \n", ++iterCount, residue);
						else ++iterCount;
					} while ((residue > accuracy) && (iterCount < maxIter));
				}
				else {
					if (verbose) std::cout << "ERROR! The cycle type should only be 0-2! \n";
					return;
				}

				solverTiming["cycles"] = cyclesTimer.get_elapsed_ms();
				solverTiming["solver_total"] = solverTotalTimer.get_elapsed_ms();
				solverTiming["iterations"] = iterCount * 1.0;
				solverTiming["residue"] = residue;
			}
			
		}
		
		else if (solverType == 4) {
			// Solve the linear system using iterative solver (Conjugate Gradient)
			Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
			cg.setTolerance(accuracy);
			cg.compute(LHS);
			{
				t0 = std::chrono::high_resolution_clock::now();
				if (x.cols() > 1)
				{
					x.col(0) = cg.solve(mx.col(0));
					t1 = std::chrono::high_resolution_clock::now();
					x.col(1) = cg.solve(mx.col(1));
					x.col(2) = cg.solve(mx.col(2));
					t2 = std::chrono::high_resolution_clock::now();
				}
				else {
					x = cg.solve(mx);
					t1 = std::chrono::high_resolution_clock::now();
				}
			}
			duration1 = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
			if (x.cols() > 1) {
				duration2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t0);
			}
		}

		Eigen::Vector3d relResidue;
		for (int i = 0; i < x.cols(); ++i) {
			Eigen::VectorXd r1 = LHS * x.col(i) - mx.col(i);
			relResidue(i) = r1.norm() / mx.col(i).norm();
		}
		double residue = relResidue.maxCoeff();
	}

	/* From SIG 2021*/
	void MultigridSolver::constructSIG21Hierarchy(const Eigen::MatrixXi& F)
	{
		std::vector<mg_data> sig21MG;
		plf::nanotimer sig21Timer;
		sig21Timer.start();
		// construct multigrid hierarchy
		int min_coarsest_nV = 500;
		float coarsening_ratio = 0.25;
		int decimation_type = 1;
		mg_precompute(V, F, coarsening_ratio, min_coarsest_nV, decimation_type, sig21MG);
		USig.resize(sig21MG.size() - 1);
		for (int k = 0; k < sig21MG.size() - 1; ++k) {
			USig[k] = sig21MG[k + 1].P;
		}
		hierarchyTiming["sig21_hierarchy"] = sig21Timer.get_elapsed_ms();
	}

	void MultigridSolver::toggleHierarchy(const Hierarchy& hierarchy)
	{
		switch (hierarchy) {
			case OURS:
				U = UOurs;
				break;
			case SIG21:
				U = USig;
				break;
			default:
				U = UOurs;
				break;
		}
	}

	void MultigridSolver::constructProlongationAblation()
	{
		// Create prolongation operator for level k+1 to k
		// Points in current level
		Eigen::MatrixXd levelPoints = V;
		Eigen::MatrixXd levelNormals = normals;
		// Points in coarser level
		Eigen::MatrixXd samplePoints;
		// Neighborhood data structure
		Eigen::MatrixXi neighLevelK = neigh;

		// Compute initial radius
		int nLevel1 = int(levelPoints.rows() / ratio);
		// if (verbose) cout << "First radius: " << radius << endl;

		// Sampling
		std::random_device					rd;
		std::default_random_engine			generator(rd());

		// Debug: set number of threads to 1
		//const int NUM_OF_THREADS = omp_get_num_procs();
		const int NUM_OF_THREADS = 1;
		omp_set_num_threads(NUM_OF_THREADS);
		int tid, ntids, ipts, istart, iproc;

		hierarchyTiming["sampling"] = 0.0;
		hierarchyTiming["cluster"] = 0.0;
		hierarchyTiming["next_neighborhood"] = 0.0;
		hierarchyTiming["next_positions"] = 0.0;
		hierarchyTiming["triangle_finding"] = 0.0;
		hierarchyTiming["triangle_selection"] = 0.0;
		hierarchyTiming["levels"] = DoF.size() - 1.0;
		// For each level
		int k = 0;
		DoF.clear();
		DoF.shrink_to_fit();
		DoF.push_back(levelPoints.rows());
		while (levelPoints.rows() > lowBound && k < 10) {
			double radius = std::cbrt(ratio) * computeAverageEdgeLength(levelPoints, neighLevelK);
			// -- Setting up variables
			// Data structure for neighbors inside level k
			std::vector<std::set<int>> neighborsList;

			// List of triplets to build prolongation operator U
			std::vector<Eigen::Triplet<double>> AllTriplet, UNeighAllTriplet;

			// The nearest coarse point for each fine point and the distances computed
			Eigen::VectorXd D(levelPoints.rows());
			D.setConstant(std::numeric_limits<double>::max());
			nearestSource.push_back(std::vector<size_t>(levelPoints.rows()));

			// -- Sample a subset for level k + 1
			if (verbose) printf("__Constructing Prolongation Operator for level = %d for ablations.\n", k);

			// Sample points that will be part of the coarser level with Poisson disk sampling
			if (verbose) std::cout << "Obtaining subset from the finer level\n";
			
			plf::nanotimer samplingTimer;
			samplingTimer.start();

			samples.push_back(fastDiskSample(levelPoints, neighLevelK, radius, D, nearestSource[k]));
			DoF.push_back(samples[k].size());

			if (samples[k].size() < lowBound) {
				nearestSource.pop_back();
				break;
			}

			if (verbose) cout << "Actual number found: " << samples[k].size() << endl;
			DoF[k + 1] = samples[k].size();

			hierarchyTiming["sampling"] += samplingTimer.get_elapsed_ms();

			plf::nanotimer clusterTimer;
			clusterTimer.start();			

			// -- Compute distance from fine points to coarse points and get closest coarse point
			// using distances from MIS if computed before
			constructDijkstraWithCluster(levelPoints, samples[k], neighLevelK, k, D, nearestSource[k]); // Stores result in nearestSource[k]
			hierarchyTiming["cluster"] += clusterTimer.get_elapsed_ms();

			plf::nanotimer nextNeighTimer;
			nextNeighTimer.start();

			// Create neighborhood for the next level
			neighborsList.resize(DoF[k + 1]);
			for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
				for (int j = 0; j < neighLevelK.cols(); ++j) {
					int neighIdx = neighLevelK(fineIdx, j);
					if (neighIdx < 0) break;
					if (nearestSource[k][fineIdx] != nearestSource[k][neighIdx]) {
						neighborsList[nearestSource[k][fineIdx]].insert(nearestSource[k][neighIdx]);
					}
				}
			}

			// Store in homogeneous datastructure
			int maxNeighNum = 0;
			for (int i = 0; i < neighborsList.size(); ++i) {
				if (neighborsList[i].size() > maxNeighNum) {
					maxNeighNum = neighborsList[i].size();
				}
			}
			neighLevelK.resize(DoF[k + 1], maxNeighNum);
			neighLevelK.setConstant(-1);
			for (int i = 0; i < neighborsList.size(); ++i) {
				neighLevelK(i, 0) = i;
				int iCounter = 1;
				for (int node : neighborsList[i]) {
					if (node == i)	continue;
					if (iCounter >= maxNeighNum) break;
					neighLevelK(i, iCounter) = node;
					iCounter++;
				}
			}
			hierarchyTiming["next_neighborhood"] += nextNeighTimer.get_elapsed_ms();

			plf::nanotimer nextPosTimer;
			nextPosTimer.start();
			if (verbose) std::cout << "Setting up the point locations for the next level\n";

			// Setting up the DoF for the next level
			// tempPoints are the centers of the voronoi cells, each row for each voronoi cells
			Eigen::MatrixXd tempPoints(DoF[k + 1], levelPoints.cols());
			tempPoints.setZero();
			if (nested) {
				for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
					tempPoints.row(coarseIdx) = levelPoints.row(samples[k][coarseIdx]);
				}
			} else {
				std::vector<int> clusterSizes(DoF[k + 1]);
				for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
					int coarseIdx = nearestSource[k][fineIdx];
					tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) + levelPoints.row(fineIdx);
					++clusterSizes[coarseIdx];
				}
				for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
					if (clusterSizes[coarseIdx] == 1) {
						tempPoints.row(coarseIdx) = levelPoints.row(samples[k][coarseIdx]);
						for (int neighIdx : neighborsList[coarseIdx]) {
							tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) + levelPoints.row(samples[k][neighIdx]);
						}
						tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) / (neighborsList[coarseIdx].size() + 1.);
					} else {
						tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) / clusterSizes[coarseIdx];
					}
				}
			}
			if (debug) levelV.push_back(tempPoints);
			hierarchyTiming["next_positions"] += nextPosTimer.get_elapsed_ms();

			plf::nanotimer triangleSelectionTimer;
			triangleSelectionTimer.start();

			// Create local point assignment
			for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
				Eigen::RowVector3d finePoint = levelPoints.row(fineIdx);
				int coarseIdx = nearestSource[k][fineIdx];
				Eigen::RowVector3d coarsePoint = tempPoints.row(coarseIdx);
				std::vector<double> weights(0);

				if (neighborsList[coarseIdx].empty()) {
					// If the coarse point has no neighbors,
					// set the weight to 1 for the coarse point.
					// Note: this should not happen.
					AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, 1.));
				} else {
					int numNeigh = neighborsList[coarseIdx].size();
					int numPoints = std::min(numNeigh, ablationNumPoints);
					std::vector<int> prolongFrom(numPoints);
					prolongFrom[0] = coarseIdx;
					if (!ablationRandom) {
						std::vector<VertexPair> pointsDistances;
						for (int j = 0; j < neighLevelK.cols(); ++j) {
							int neighIdx = neighLevelK(coarseIdx, j);
							if (neighIdx < 0 || neighIdx == coarseIdx) continue;
							VertexPair vp = {neighIdx, (finePoint - tempPoints.row(neighIdx)).norm()};
							pointsDistances.push_back(vp);
						}
						std::sort(pointsDistances.begin(), pointsDistances.end(), std::less<VertexPair>());
						for (int j = 1; j < numPoints; ++j) {
							prolongFrom[j] = pointsDistances[j - 1].vId;
						}
					} else {
						std::vector<int> randomNeigh(numNeigh);
						std::iota(randomNeigh.begin(), randomNeigh.end(), 0);
						std::shuffle(randomNeigh.begin(), randomNeigh.end(), generator);
						for (int j = 1; j < numPoints; ++j) {
							prolongFrom[j] = neighLevelK(coarseIdx, randomNeigh[j - 1]);
						}
					}
					weights = inverseDistanceWeights(tempPoints, finePoint, prolongFrom);
					for (int j = 0; j < prolongFrom.size(); j++) {
						AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, prolongFrom[j], weights[j]));
					}
				}
			}
			hierarchyTiming["triangle_selection"] += triangleSelectionTimer.get_elapsed_ms();

			levelPoints = tempPoints;

			Eigen::SparseMatrix<double> ULevel;
			ULevel.resize(DoF[k], DoF[k + 1]);
			ULevel.setFromTriplets(AllTriplet.begin(), AllTriplet.end());
			U.push_back(ULevel);
			AllTriplet.clear();
			AllTriplet.shrink_to_fit();
			++k;
		}
	}
}

