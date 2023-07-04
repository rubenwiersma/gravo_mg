#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "plf_nanotimer.h"

#include <set>
#include <vector>
#include <random>

#include <Eigen/Eigen>

using namespace std;

namespace MGBS {
	void scaleMesh(Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double scaleRatio = 1.0);
	void normalize_unit_area(Eigen::MatrixXd& V,const Eigen::MatrixXi& F);
	void writeSparseMatrixToFile(const Eigen::SparseMatrix<double>& M, const string& filename);
	void writeMatrixToFile(const Eigen::MatrixXd& M, const string& filename);
	void textureFromPNG(const Eigen::VectorXd& f, Eigen::MatrixXd& FColor);
	void writeTiming(const map<string, double>& timing, const string& experiment, const string& filename, const bool& writeHeaders);
	void writeConvergence(const std::vector<std::tuple<double, double>>& convergence, const string& filename);
}

#endif // !UTILITY_H