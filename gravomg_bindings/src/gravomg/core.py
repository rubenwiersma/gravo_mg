import numpy as np
from  scipy.sparse import csr_matrix
import gravomg_bindings

from gravomg_bindings import Hierarchy, Sampling, Weighting

class MultigridSolver(object):
    def __init__(
        self, pos, neigh, mass,
        ratio=8.0, lower_bound=1000, cycle_type=0, tolerance=1e-4, stopping_criteria=2, pre_iters=2, post_iters=2, max_iter=100, 
        check_voronoi=True, nested=False, sampling_strategy=Sampling.FASTDISK, weighting=Weighting.BARYCENTRIC,
        sig06=False, normals=None, verbose=False, debug=False, ablation=False, ablation_num_points=3, ablation_random=False
        ):
        """Creates the Gravo MG solver for linear systems on curved surfaces (meshes and point clouds).
        
        Args:
            pos (np.ndarray): The positions of the points in the mesh or point cloud.
            neigh (np.ndarray): The neighbors of each point in the mesh or point cloud.
                This should be given as a homogeneous array of size (n_points, max_neighbors),
                padded with -1.
            mass (scipy.sparse.csr_matrix): The mass matrix of the mesh or point cloud.
            ratio (float, optional): The coarsening ratio. Defaults to 8.0.
            lower_bound (int, optional): The lower bound on the number of points in the coarsest level.
                Defaults to 1000.
            cycle_type (int, optional): The type of cycle to use. Defaults to 0 (V-cycle).
            tolerance (float, optional): The tolerance for the stopping criteria. Defaults to 1e-4.
            stopping_criteria (int, optional): The stopping criteria to use. Defaults to 2 (relative residual).
            pre_iters (int, optional): The number of pre-smoothing iterations. Defaults to 2.
            post_iters (int, optional): The number of post-smoothing iterations. Defaults to 2.
            max_iter (int, optional): The maximum number of iterations. Defaults to 100.
            check_voronoi (bool, optional): Whether to use the Voronoi diagram to compose candidate triangles.
                Defaults to True.
            nested (bool, optional): When set to True, does not shift the coarser points to the barycenter.
                Defaults to False.
            sampling_strategy (int, optional): The sampling strategy to use. Defaults to Sampling.FASTDISK.
            weighting (int, optional): The weighting scheme to use. Defaults to Weighting.BARYCENTRIC.
            sig06 (bool, optional): Whether to use the SIG06 hierarchy. Defaults to False.
            normals (np.ndarray, optional): The normals of the points in the mesh or point cloud.
                Defaults to None.
            verbose (bool, optional): Whether to print verbose output. Defaults to False.
            debug (bool, optional): Whether to print debug output. Defaults to False.
            ablation (bool, optional): Whether to use the ablation hierarchy. Defaults to False.
            ablation_num_points (int, optional): The number of points to use in the ablation hierarchy.
                Defaults to 3.
            ablation_random (bool, optional): Whether to use random points in the ablation hierarchy.
                Defaults to False.
        """
        super().__init__()
        if not mass.getformat() == 'csr':
            mass = mass.tocsr()
        normals = pos if normals is None else normals
        self.solver = gravomg_bindings.MultigridSolver(
            pos, neigh, mass,
            ratio, lower_bound, cycle_type, tolerance, stopping_criteria, pre_iters, post_iters, max_iter,
            check_voronoi, nested, sampling_strategy, weighting,
            sig06, normals, verbose, debug, ablation, ablation_num_points, ablation_random
            )
        self.sig21_computed = False
        self.sig21bary_computed = False

    def construct_sig21_hierarchy(self, faces):
        """Constructs the hierarchy as proposed by Liu et al. [2021]
        with the default parameters given in their paper.

        Args:
            faces (np.ndarray): The faces of the mesh.
        """
        self.sig21_computed = True
        self.solver.construct_sig21_hierarchy(faces)

    def toggle_hierarchy(self, hierarchy_type):
        """Toggles the hierarchy to use between our method and Liu et al. [2021]."""
        assert (
            hierarchy_type == Hierarchy.OURS
            or (hierarchy_type == Hierarchy.SIG21 and self.sig21_computed)
            or (hierarchy_type == Hierarchy.SIG21BARY and self.sig21bary_computed)
        )
        self.solver.toggle_hierarchy(hierarchy_type)

    def solve(self, lhs, rhs):
        """Solves a linear system Ax = b, where lhs is A and rhs is b.

        Args:
            lhs (scipy.sparse.csr_matrix): The left-hand side of the linear system.
            rhs (np.ndarray): The right-hand side of the linear system.
        """
        if not lhs.getformat() == 'csr':
            print('LHS is not in CSR format, converting to CSR')
            lhs = lhs.tocsr()
        return self.solver.solve(lhs, rhs)

    def direct_solve(self, lhs, rhs, pardiso=False):
        """Solves a linear system with a direct solver.
        If pardiso is set to True, uses the Pardiso solver.
        """
        return self.solver.direct_solve(lhs, rhs, pardiso)

    # Getters and setters

    @property
    def prolongation_matrices(self):
        return self.solver.prolongation_matrices()

    def set_prolongation_matrices(self, U):
        self.solver.set_prolongation_matrices(U)

    @property
    def sampling_indices(self):
        return self.solver.sampling_indices()

    @property
    def level_points(self):
        return self.solver.level_points()

    @property
    def level_edges(self):
        return self.solver.level_edges()

    @property
    def notrimap(self):
        return self.solver.notrimap()

    @property
    def all_triangles(self):
        return self.solver.all_triangles()

    @property
    def coarse_normals(self):
        return self.solver.coarse_normals()

    @property
    def nearest_source(self):
        return self.solver.nearest_source()

    # Functions to write timing logs to a file.
    
    def write_hierarchy_timing(self, experiment, file, write_headers=False):
        return self.solver.write_hierarchy_timing(experiment, file, write_headers)

    def write_solver_timing(self, experiment, file, write_headers=False):
        return self.solver.write_solver_timing(experiment, file, write_headers)

    def write_convergence(self, file):
        return self.solver.write_convergence(file)

    def residual(self, lhs, rhs, solution, type=2):
        return self.solver.residual(lhs, rhs, solution, type)