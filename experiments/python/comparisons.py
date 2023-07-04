import pathlib, argparse
import time

import igl
from robust_laplacian import mesh_laplacian, point_cloud_laplacian

import numpy as np
from scipy import sparse

from gravomg import MultigridSolver, Hierarchy, Sampling, Weighting
from gravomg.util import neighbors_from_stiffness, normalize_area, normalize_axes, normalize_bounding_box

from util import read_mesh, read_pointcloud
from comparisons_to_table import save_to_table

from pyamg import ruge_stuben_solver, smoothed_aggregation_solver

pyamg_iterations = 0

def list_shapes(dir_path):
    shape_files = []

    d = pathlib.Path(dir_path)
    for entry in d.iterdir():
        if entry.is_file() and (entry.match('*.obj') or entry.match('*.ply') or entry.match('*.off')):
            shape_files.append(entry)

    return shape_files

def preprocess(V, args, F=None):
    # Normalize shape
    if not args.pointcloud:
        V = normalize_area(V, F)
        N = igl.per_vertex_normals(V, F)
    else:
        V = normalize_bounding_box(V)
        N = None

    # Compute operators
    if not args.robust and not args.pointcloud:
        S = -igl.cotmatrix(V, F)
        M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_VORONOI)
    elif args.pointcloud:
        S, M = point_cloud_laplacian(V)
    elif args.robust:
        S, M = mesh_laplacian(V, F)

    Minv = sparse.diags(1 / M.diagonal())
    if args.robust_neigh:
        S_robust, _ = mesh_laplacian(V, F)
        neigh = neighbors_from_stiffness(S_robust)
    else:
        neigh = neighbors_from_stiffness(S)
    B = S @ Minv @ S
    return V, F, N, M, S, neigh, B

def smoothe_dir(args):
    global pyamg_iterations

    shape_files = list_shapes(args.in_dir)
    print(f'{len(shape_files)} files found in folder \'{args.in_dir}\'')

    print('Starting experiments...')
    for i, shape_file in enumerate(shape_files):
        print(f'Shape {i + 1}/{len(shape_files)}: {shape_file.stem}')
        if not args.pointcloud:
            V, F = read_mesh(shape_file)
        else:
            V = read_pointcloud(shape_file)
            F = None

        V, F, N, M, S, neigh, B = preprocess(V, args, F)
       
        # Setup system
        if args.poisson:
            lhs = M * args.tau + S if not args.bilaplacian else M * args.tau + B
        else:
            lhs = M + args.tau * S if not args.bilaplacian else M + args.tau * B
        lhs = lhs.tocsr()

        # Create our hierarchy to  use points for SIG21
        solver = MultigridSolver(V, neigh, M, ratio=args.ratio, lower_bound=args.lower_bound, tolerance=args.tolerance, nested=args.nested, sampling_strategy=args.sampling, verbose=args.verbose)

        # Right hand side: random
        rng = np.random.default_rng(seed=args.seed)
        if args.input_smooth:
            max_idx = np.argmax(V.sum(axis=1))
            min_idx = np.argmin(V.sum(axis=1))
            y = np.zeros_like(V[:, 0:1])
            y[max_idx] = 1
            y[min_idx] = -1
            y = solver.solve(M + 0.5 * S, M @ y).flatten()
            y = y + rng.standard_normal(V.shape[0]) * 0.0000005
        else:
            y = rng.standard_normal((V.shape[0], 1))
        rhs = M @ y

        if (args.direct):
            print('Direct solver')
            solver.direct_solve(lhs, rhs)
            solver.direct_solve(lhs, rhs, pardiso=True)
            solver.write_solver_timing(shape_file.stem, f'{args.out_dir}/direct_tau{args.tau}_{args.label}.csv', write_headers=i==0)

        # Only run SIG21 code for meshes it can work on
        if (args.sig21):
            print('Constructing and solving SIG21')
            solver.construct_sig21_hierarchy(F)
            solver.write_hierarchy_timing(shape_file.stem, f'{args.out_dir}/hierarchy_sig21_{args.label}.csv', write_headers=i==0)

            # Solve SIGGRAPH21
            solver.toggle_hierarchy(Hierarchy.SIG21)
            solver.solve(lhs, rhs)
            solver.write_solver_timing(shape_file.stem, f'{args.out_dir}/solver_sig21_tau{args.tau}_{args.label}.csv', write_headers=i==0)
            solver.write_convergence(f'{args.out_dir}/convergence/sig21/{shape_file.stem}_tau{args.tau}_{args.label}.csv')

        if (args.sig06):
            print('Constructing and solving SIG06')
            sig06_solver = MultigridSolver(V, neigh, M, sig06=True, ratio=args.ratio, lower_bound=args.lower_bound, tolerance=args.tolerance)
            sig06_solver.write_hierarchy_timing(shape_file.stem, f'{args.out_dir}/hierarchy_sig06_{args.label}.csv', write_headers=i==0)

            # Solve SIGGRAPH06
            sig06_solver.solve(lhs, rhs)
            sig06_solver.write_solver_timing(shape_file.stem, f'{args.out_dir}/solver_sig06_tau{args.tau}_{args.label}.csv', write_headers=i==0)
            sig06_solver.write_convergence(f'{args.out_dir}/convergence/sig06/{shape_file.stem}_tau{args.tau}_{args.label}.csv')

        if (args.amg):
            pyamg_iterations = 0
            def count_iterations(xk):
                global pyamg_iterations
                residual = solver.residual(lhs, rhs[:, 0], xk)
                if residual > args.tolerance:
                    pyamg_iterations += 1
            print('Constructing and solving PyAMG: Ruge-Stuben')
            t = time.perf_counter()
            rs_solver = ruge_stuben_solver(lhs)
            rs_hierarchy_time = time.perf_counter() - t

            rs_solver.solve(rhs[:, 0], tol=1e-12, callback=count_iterations) 
            t = time.perf_counter()
            rs_solver.solve(rhs[:, 0], tol=1e-12, maxiter=pyamg_iterations) 
            rs_solve_time = time.perf_counter() - t

            # Write timings to CSV
            rs_file = pathlib.Path(args.out_dir) / f'amg_rs_tau{args.tau}_{args.label}.csv'
            with rs_file.open('w' if i == 0 else 'a') as f:
                if i == 0: f.write('experiment,rs_hierarchy,rs_iterations,rs_solver\n')
                f.write(f'{shape_file.stem},{rs_hierarchy_time},{pyamg_iterations},{rs_solve_time}\n')

            pyamg_iterations = 0
            print('Constructing and solving PyAMG: Smoothed Aggregation')
            t = time.perf_counter()
            sa_solver = smoothed_aggregation_solver(lhs)
            sa_hierarchy_time = time.perf_counter() - t

            sa_solver.solve(rhs[:, 0], tol=1e-8, callback=count_iterations) 
            t = time.perf_counter()
            sa_solver.solve(rhs[:, 0], tol=1e-8, maxiter=pyamg_iterations) 
            sa_solve_time = time.perf_counter() - t

            # Write timings to CSV
            sa_file = pathlib.Path(args.out_dir) / f'amg_sa_tau{args.tau}_{args.label}.csv'
            with sa_file.open('w' if i == 0 else 'a') as f:
                if i == 0: f.write('experiment,sa_hierarchy,sa_iterations,sa_solver\n')
                f.write(f'{shape_file.stem},{sa_hierarchy_time},{pyamg_iterations},{sa_solve_time}\n')


        for j in range(args.num_repetitions): 
            print(f'Constructing and solving OURS ({j + 1}/{args.num_repetitions})')
            solver = MultigridSolver(V, neigh, M, normals=N, ratio=args.ratio, lower_bound=args.lower_bound, check_voronoi=not args.all_triangles, tolerance=args.tolerance, nested=args.nested, sampling_strategy=args.sampling, weighting=args.weighting, ablation=args.ablation, ablation_num_points=args.ablation_n, ablation_random=args.ablation_random)
            solver.write_hierarchy_timing(shape_file.stem, f'{args.out_dir}/hierarchy_ours_{args.label}.csv', write_headers=(i==0 and j==0))
            solver.toggle_hierarchy(Hierarchy.OURS)
            solver.solve(lhs, rhs)
            solver.write_solver_timing(shape_file.stem, f'{args.out_dir}/solver_ours_tau{args.tau}_{args.label}.csv', write_headers=(i==0 and j==0))
            solver.write_convergence(f'{args.out_dir}/convergence/ours/{shape_file.stem}_tau{args.tau}_{args.label}.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run MultigridSolver benchmark on directory')
    parser.add_argument('--tau', type=float, default=1e-3, help='Smoothing coefficient (default: 1e-3)')
    parser.add_argument('--ratio', type=float, default=8, help='Decay of density for hierarchy construction (default: 2)')
    parser.add_argument('--lower_bound', type=int, default=1000, help='Minimum number of points in the lowest level (default: 500)')
    parser.add_argument('--tolerance', type=float, default=1e-4, help='Desired solving accuracy (default: 1e-4)')
    parser.add_argument('--label', type=str, default='laplacian', help='Label to use for output data (default: laplacian)')
    parser.add_argument('--in_dir', type=str, default='../../data/robust', help='Directory to run benchmark on, relative to this script (default: ../../data/final)')
    parser.add_argument('--out_dir', type=str, default='../../out/timing', help='Directory to output data, relative to this script (default: ../../out/timing)')
    parser.add_argument('--num_repetitions', type=int, default=1, help='Number of repeated constructions and solves to averatge (default: 11)')
    parser.add_argument('--bilaplacian', action='store_true', help='Solve a Bilaplacian system')
    parser.add_argument('--poisson', action='store_true', help='Solve a Poisson problem')
    parser.add_argument('--input_smooth', action='store_true', help='Uses a relatively smooth input function with little noise added')
    parser.add_argument('--robust', action='store_true', help='Uses robust Laplacian if provided')
    parser.add_argument('--robust_neigh', action='store_true', help='Uses robust Laplacian neighborhoods, but cotan Laplacian for linear system if provided')
    parser.add_argument('--pointcloud', action='store_true', help='Uses only point positions if provided')
    parser.add_argument('--nonmanifold', action='store_true', help='Uses robust Laplacian and only runs our method')
    parser.add_argument('--all_triangles', action='store_true', help='Considers all triangles, instead of only Voronoi neighbors')
    parser.add_argument('--nested', action='store_true', help='Maintains point positions throughout hierarchy')
    parser.add_argument('--direct', action='store_true', help='Includes direct solver')
    parser.add_argument('--nosig21', action='store_true', help='Includes SIG21 results')
    parser.add_argument('--sig06', action='store_true', help='Includes SIG06 results')
    parser.add_argument('--amg', action='store_true', help='Includes PyAMG results')
    parser.add_argument('--sampling', type=str, default='fastdisk', choices=['fastdisk', 'poissondisk', 'random', 'fps', 'mis'], help='Sampling method to use (default: poissondisk)')
    parser.add_argument('--weighting', type=str, default='barycentric', choices=['barycentric', 'uniform', 'invdist'], help='Sampling method to use (default: poissondisk)')
    parser.add_argument('--ablation', action='store_true', help='Picks closest n points')
    parser.add_argument('--ablation_n', type=int, default=3, help='The number of points to select')
    parser.add_argument('--ablation_random', action='store_true', help='Picks random points in the voronoi neighborhood')
    parser.add_argument('--no_names', action='store_true', help='Removes first two columns from latex output')
    parser.add_argument('--verbose', action='store_true', help='Prints debug information to the console')
    parser.add_argument('--seed', type=int, default=42, help='Seed to use for random function generation')
    args = parser.parse_args()

    args.in_dir = str((pathlib.Path(__file__).parent.absolute() / args.in_dir).resolve())
    args.out_dir = str((pathlib.Path(__file__).parent.absolute() / args.out_dir).resolve())

    sampling_enums = {'fastdisk': Sampling.FASTDISK, 'poissondisk': Sampling.POISSONDISK, 'random': Sampling.RANDOM, 'fps': Sampling.FPS, 'mis': Sampling.MIS}
    args.sampling = sampling_enums[args.sampling]
    weighting_enums = {'barycentric': Weighting.BARYCENTRIC, 'uniform': Weighting.UNIFORM, 'invdist': Weighting.INVDIST}
    args.weighting = weighting_enums[args.weighting]
    args.sig21 = not (args.pointcloud or args.nonmanifold) and not args.nosig21

    experiment_details = args.label + '\n--\nSettings:\n--\n'
    for arg in vars(args):
        experiment_details += '{}: {}\n'.format(arg, getattr(args, arg))
    with (pathlib.Path(args.out_dir) / f'settings_{args.label}_tau{args.tau}.txt').open('w') as f:
        f.write(experiment_details)

    print(experiment_details)
    print('---')

    smoothe_dir(args)
    save_to_table(args.out_dir, args.tau, args.label, sig21=args.sig21, sig06=args.sig06, amg=args.amg, direct=args.direct, names_counts=not args.no_names)