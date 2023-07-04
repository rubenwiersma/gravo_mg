#include "gravomg/utility.h"

#include <igl/doublearea.h>

#include <fstream>


namespace MGBS {
	void scaleMesh(Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double scaleRatio)
	{
		//cout << "> Scaling the mesh to have a unit length \n";
		Eigen::Vector3d minV;
		Eigen::Vector3d maxV;
		Eigen::Vector3d length;
		Eigen::MatrixXd MatSubs(V.rows(), V.cols());
		Eigen::MatrixXd MatAdd(V.rows(), V.cols());
		double maxVal;
		double scalingRatio = scaleRatio; //dia-meter, not radius


		/* Get the min and max coefficients */
		for (int i = 0; i < V.cols(); i++)
		{
			minV(i) = V.col(i).minCoeff();
			maxV(i) = V.col(i).maxCoeff();
			length(i) = maxV(i) - minV(i);
			MatSubs.col(i).setConstant(minV(i));
		}

		maxVal = length.maxCoeff();

		/* Translate to the Origin */
		V = V - MatSubs;

		/* Scale w.r.t the longest axis */
		V = V * (scalingRatio / maxVal);

		for (int i = 0; i < V.cols(); i++) {
			maxV(i) = V.col(i).maxCoeff();
			MatAdd.col(i).setConstant(0.5 * maxV(i));
		}

		/* Translate s.t. the center is in the Origin */
		V = V - MatAdd;

		for (int i = 0; i < V.cols(); i++)
		{
			minV(i) = V.col(i).minCoeff();
			maxV(i) = V.col(i).maxCoeff();
			length(i) = maxV(i) - minV(i);
		}
		maxVal = length.maxCoeff();
	}

	void normalize_unit_area(Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
	{
		Eigen::VectorXd FA;
		igl::doublearea(V, F, FA);

		double scale = sqrt(FA.sum() / 2);
		V /= scale;
				
		V.col(0) = V.col(0).array() - V.col(0).mean();
		V.col(1) = V.col(1).array() - V.col(1).mean();
		//V.col(2) = V.col(2).array() - V.col(2).minCoeff();
		V.col(2) = V.col(2).array() - V.col(2).mean();

		std::cout << "Check mass matrix after normalization \n";
		igl::doublearea(V, F, FA);
		std::cout << "Area: " << FA.sum() / 2 << " and M(0,0)=" << FA(0) / 2 << std::endl;
	}

	void writeSparseMatrixToFile(const Eigen::SparseMatrix<double>& M, const string& filename)
	{
		ofstream myfile(filename.c_str());
		if (myfile.is_open())
		{
			cout << "Write file to text \n";
			for (int k = 0; k < M.outerSize(); ++k) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {					
					myfile << it.row() << "\t" << it.col() << "\t" << it.value() << "\n";
				}
			}
			myfile.close();
		}
		else cout << "Unable to open file: " << filename << endl;
	}

	void writeMatrixToFile(const Eigen::MatrixXd& M, const string& filename)
	{
		ofstream myfile(filename.c_str());
		if (myfile.is_open())
		{
			cout << "Write file to text \n";
			for (int i = 0; i < M.rows(); ++i) {
				for (int j = 0; j < M.cols(); ++j) {					
					myfile << M(i, j) << "\t"; 
				}
				myfile << "\n";
			}
			myfile.close();
		}
		else cout << "Unable to open file: " << filename << endl;
	}

	void writeTiming(const map<string, double>& timing, const string& experiment, const string& filename, const bool& writeHeaders = false)
	{
		ofstream myfile;
		if (writeHeaders) 
			myfile.open(filename.c_str());
		else
			myfile.open(filename.c_str(), ios_base::app);

		if (myfile.is_open())
		{
			if (writeHeaders) {
				myfile << "experiment";
				for (auto const& t : timing) {
					myfile << ',' << t.first;
				}
				myfile << "\n";
			}
			myfile << experiment;
			for (auto const& t : timing) {
				myfile << ',' << t.second;
			}
			myfile << "\n";
			myfile.close();
		}
		else cout << "Unable to open timing file: " << filename << endl;
	}	

	void writeConvergence(const std::vector<std::tuple<double, double>>& convergence, const string& filename)
	{
		ofstream myfile;
		myfile.open(filename.c_str());

		if (myfile.is_open())
		{
			myfile << "time,residue\n";
			for (auto const& it : convergence) {
				double time, residue;
				std::tie(time, residue) = it;
				myfile << time << ',' << residue << "\n";
			}
			myfile.close();
		}
		else cout << "Unable to open convergence file: " << filename << endl;
	}	

}
