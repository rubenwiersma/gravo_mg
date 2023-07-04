#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <Eigen/Dense>

/* Data structure for Priority Queue */
struct VertexPair {
	int		vId;
	double	distance;
	bool	operator> (const VertexPair& ref) const { return distance > ref.distance; }
	bool	operator< (const VertexPair& ref) const { return distance < ref.distance; }
};

std::vector<int> constructFarthestPointSample(const Eigen::MatrixXd& points,
											 const size_t numSamples,
											 const Eigen::MatrixXi& neigh);

void computeDijkstra(const Eigen::MatrixXd& points, int source, const Eigen::MatrixXi& neigh, Eigen::VectorXd &D);

#endif // !SAMPLING_H
