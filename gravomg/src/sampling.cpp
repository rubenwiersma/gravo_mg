#include "gravomg/sampling.h"

#include <queue>
#include <random>

std::vector<int> constructFarthestPointSample(const Eigen::MatrixXd& points,
											 const size_t numSamples,
											 const Eigen::MatrixXi& neigh)
{
	/* */
	std::vector<int> sampleID(numSamples);

	Eigen::VectorXd D(points.rows());
	D.setConstant(std::numeric_limits<double>::infinity());

	std::random_device rd;									// Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd());									// Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dist(0, points.rows() - 1);	// From 0 to (number of points - 1)
	int seed = dist(gen);

	sampleID[0] = seed;

	for (size_t i = 1; i < numSamples; i++) {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstra(points, sampleID[i-1], neigh, D);
		D.maxCoeff(&maxIndex);
		sampleID[i] = maxIndex;
	}

	return sampleID;
}

// Performing Dijkstra computation on point clouds
void computeDijkstra(const Eigen::MatrixXd& points, int source, const Eigen::MatrixXi& neigh, Eigen::VectorXd &D)
{
	std::priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistanceQueue;

	D(source) = 0.0;
	VertexPair vp{ source, D(source) };
	DistanceQueue.push(vp);

	//double distanceFromSource = std::numeric_limits<double>::infinity();

	while (!DistanceQueue.empty()){
		VertexPair vp1 = DistanceQueue.top();
		Eigen::RowVector3d vertex1 = points.row(vp1.vId);
		DistanceQueue.pop();

		for (int i = 0; i < neigh.cols(); ++i){
			int vNeigh = neigh(vp1.vId, i);

			if (vNeigh >= 0) {
				double dist, distTemp;
				Eigen::RowVector3d vertex2 = points.row(vNeigh);
				dist = (vertex2 - vertex1).norm();
				distTemp = vp1.distance + dist;
				if(distTemp < D(vNeigh)){
					D(vNeigh) = distTemp;
					VertexPair v2{ vNeigh, distTemp };
					DistanceQueue.push(v2);
				}
			}
		}
	}

}
