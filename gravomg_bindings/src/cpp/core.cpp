#include "gravomg/multigrid_solver.h"
#include "gravomg/utility.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/Sparse"

namespace py = pybind11;

class MultigridSolver {

public:
    // Expects positions as N x 3 matrix and neighbors as N x K,
    // where K is the maximum number of neighbors.
    // Row i contains the indices of the neighbors of node i;
    // the neighbors should be padded with -1 to get K entries per row.
    MultigridSolver(
        Eigen::MatrixXd positions, Eigen::MatrixXi neighbors, Eigen::SparseMatrix<double> mass,
        double ratio, int low_bound, // Hierarchy settings
        int cycle_type, double tolerance, int stopping_criteria, int pre_iters, int post_iters, int max_iter, // Solver settings
        bool check_voronoi, bool nested, Sampling sampling_strategy, Weighting weighting, bool sig06, Eigen::MatrixXd normals, bool verbose, bool debug, // Debug settings
        bool ablation, int ablation_num_points, bool ablation_random // Ablations
        ) {
        solver.reset(new MGBS::MultigridSolver(positions, neighbors, mass));
        solver->checkVoronoi = check_voronoi;
        solver->nested = nested;
        solver->samplingStrategy = sampling_strategy;
        solver->weightingScheme = weighting;
        solver->maxIter = max_iter;
        solver->sig06 = sig06;

        // Ablation for finding points
        solver->ablation = ablation;
        solver->ablationNumPoints = ablation_num_points;
        solver->ablationRandom = ablation_random;

        solver->normals = normals;

        solver->verbose = verbose;
        solver->debug = debug;

        // Building hierarchy
        // solver->setSamplingNumber(1, -1, ratio * ratio, true, low_bound);
        solver->ratio = ratio;
        solver->lowBound = low_bound;
        solver->buildHierarchy();

        // Set solver settings
        solver->cycleType = cycle_type;
        solver->accuracy = tolerance;
        solver->stoppingCriteria = stopping_criteria;
        solver->preIters = pre_iters;
        solver->postIters = post_iters;
        solver->isSmootherGaussSeidel = true;
    }

    void construct_sig21_hierarchy(Eigen::MatrixXi F) {
        solver->constructSIG21Hierarchy(F);
    }

    void toggle_hierarchy(Hierarchy hierarchy) {
        solver->toggleHierarchy(hierarchy);
    }

    Eigen::MatrixXd solve(Eigen::SparseMatrix<double> lhs, Eigen::MatrixXd rhs) {
        Eigen::MatrixXd x = rhs;
        solver->solve(lhs, rhs, x, 2);
        return x;
    }

    Eigen::MatrixXd direct_solve(Eigen::SparseMatrix<double> lhs, Eigen::MatrixXd rhs, bool pardiso) {
        Eigen::MatrixXd x = rhs;
        solver->solve(lhs, rhs, x, (pardiso) ? 1 : 0);
        return x;
    }

    //-- Data access

    std::vector<Eigen::SparseMatrix<double>> prolongation_matrices() {
        return solver->U;
    }

    void set_prolongation_matrices(std::vector<Eigen::SparseMatrix<double>> U) {
        solver->U = U;
    }

    std::vector<std::vector<int>> sampling_indices() {
        return solver->samples;
    }

    std::vector<Eigen::MatrixXd> level_points() {
        return solver->levelV;
    }
    
    std::vector<Eigen::MatrixXi> level_edges() {
        return solver->levelE;
    }

    std::vector<std::vector<int>> notrimap() {
        return solver->noTriFoundMap;
    }

    std::vector<std::vector<std::vector<int>>> all_triangles() {
        return solver->allTriangles;
    }

    std::vector<Eigen::MatrixXd> coarse_normals() {
        return solver->levelN;
    }

    std::vector<std::vector<size_t>> nearest_source() {
        return solver->nearestSource;
    }

    void write_hierarchy_timing(string experiment, string file, bool write_headers) {
        MGBS::writeTiming(solver->hierarchyTiming, experiment, file, write_headers);
    }

    void write_solver_timing(string experiment, string file, bool write_headers) {
        MGBS::writeTiming(solver->solverTiming, experiment, file, write_headers);
    }

    void write_convergence(string file) {
        MGBS::writeConvergence(solver->convergence, file);
    }

    // -- Metrics

    double residual(Eigen::SparseMatrix<double> lhs, Eigen::MatrixXd rhs, Eigen::MatrixXd solution, int type=2) {
        return solver->residualCheck(lhs, rhs, solution, type);
    }


private:
    std::unique_ptr<MGBS::MultigridSolver> solver; 
};


PYBIND11_MODULE(gravomg_bindings, m) {
    m.doc() = "Multigrid solver bindings";

    py::class_<MultigridSolver>(m, "MultigridSolver")
        .def(py::init<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::SparseMatrix<double>, double, int, int, double, int, int, int, int, bool, bool, Sampling, Weighting, bool, Eigen::MatrixXd, bool, bool, bool, int, bool>())
        .def("construct_sig21_hierarchy", &MultigridSolver::construct_sig21_hierarchy, py::arg("F"))
        .def("toggle_hierarchy", &MultigridSolver::toggle_hierarchy, py::arg("hierarchy"))
        .def("solve", &MultigridSolver::solve, py::arg("lhs"), py::arg("rhs"))
        .def("direct_solve", &MultigridSolver::direct_solve, py::arg("lhs"), py::arg("rhs"), py::arg("pardiso"))
        .def("prolongation_matrices", &MultigridSolver::prolongation_matrices)
        .def("set_prolongation_matrices", &MultigridSolver::set_prolongation_matrices, py::arg("U"))
        .def("sampling_indices", &MultigridSolver::sampling_indices)
        .def("level_points", &MultigridSolver::level_points)
        .def("level_edges", &MultigridSolver::level_edges)
        .def("notrimap", &MultigridSolver::notrimap)
        .def("all_triangles", &MultigridSolver::all_triangles)
        .def("coarse_normals", &MultigridSolver::coarse_normals)
        .def("nearest_source", &MultigridSolver::nearest_source)
        .def("write_hierarchy_timing", &MultigridSolver::write_hierarchy_timing, py::arg("experiment"), py::arg("file"), py::arg("write_headers"))
        .def("write_solver_timing", &MultigridSolver::write_solver_timing, py::arg("experiment"), py::arg("file"), py::arg("write_headers"))
        .def("write_convergence", &MultigridSolver::write_convergence, py::arg("file"))
        .def("residual", &MultigridSolver::residual, py::arg("lhs"), py::arg("rhs"), py::arg("solution"), py::arg("type") = 2);

    py::enum_<Hierarchy>(m, "Hierarchy")
        .value("OURS", OURS)
        .value("SIG21", SIG21);
        
    py::enum_<Sampling>(m, "Sampling")
        .value("FASTDISK", FASTDISK)
        .value("POISSONDISK", POISSONDISK)
        .value("FPS", FPS)
        .value("RANDOM", RANDOM)
        .value("MIS", MIS);

    py::enum_<Weighting>(m, "Weighting")
        .value("BARYCENTRIC", BARYCENTRIC)
        .value("UNIFORM", UNIFORM)
        .value("INVDIST", INVDIST);
}
