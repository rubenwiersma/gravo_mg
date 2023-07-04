import pandas as pd
from tabulate import tabulate

from pathlib import Path

def float_to_string(float_col, float_format='{:.2f}'):
    return float_col.map(float_format.format)

def save_to_table(out_dir, tau, label, latex=True, sig21=True, sig06=False, amg=False, direct=False, std=False, names_counts=True):
    hierarchy_path = f'{out_dir}/hierarchy_ours_{label}.csv'
    solver_ours_path = f'{out_dir}/solver_ours_tau{tau}_{label}.csv'
    hierarchy_data = pd.read_csv(hierarchy_path).sort_values('experiment').reset_index()
    solver_ours_data = pd.read_csv(solver_ours_path).sort_values('experiment').reset_index()

    if direct:
        solver_direct_path = f'{out_dir}/direct_tau{tau}_{label}.csv'
        solver_direct_data = pd.read_csv(solver_direct_path).sort_values('experiment').reset_index()

    if sig21:
        hierarchy_sig21_path = f'{out_dir}/hierarchy_sig21_{label}.csv'
        solver_sig21_path = f'{out_dir}/solver_sig21_tau{tau}_{label}.csv'
        hierarchy_sig21_data = pd.read_csv(hierarchy_sig21_path).sort_values('experiment').reset_index()
        solver_sig21_data = pd.read_csv(solver_sig21_path).sort_values('experiment').reset_index().rename(
            columns={'iterations': 'sig21_iterations', 'residue': 'sig21_residue', 'solver_total': 'sig21_solver'}
        )

    if sig06:
        hierarchy_sig06_path = f'{out_dir}/hierarchy_sig06_{label}.csv'
        solver_sig06_path = f'{out_dir}/solver_sig06_tau{tau}_{label}.csv'
        hierarchy_sig06_data = pd.read_csv(hierarchy_sig06_path).sort_values('experiment').reset_index().rename(
            columns={'hierarchy': 'sig06_hierarchy'}
        )
        solver_sig06_data = pd.read_csv(solver_sig06_path).sort_values('experiment').reset_index().rename(
            columns={'iterations': 'sig06_iterations', 'residue': 'sig06_residue', 'solver_total': 'sig06_solver'}
        )

    if amg:
        amg_rs_path = f'{out_dir}/amg_rs_tau{tau}_{label}.csv'
        amg_rs_data = pd.read_csv(amg_rs_path).sort_values('experiment').reset_index()
        amg_sa_path = f'{out_dir}/amg_sa_tau{tau}_{label}.csv'
        amg_sa_data = pd.read_csv(amg_sa_path).sort_values('experiment').reset_index()

    hierarchy_data_summaries = hierarchy_data.groupby('experiment', as_index=False).agg(
        n_vertices = ('n_vertices', 'max'),
        mean_hierarchy = ('hierarchy', 'mean'),
        std_hierarchy = ('hierarchy', 'std')
    )

    solver_ours_summaries = solver_ours_data.groupby('experiment', as_index=False).agg(
        median_iterations = ('iterations', 'median'),
        mean_iterations = ('iterations', 'mean'),
        std_iterations = ('iterations', 'std'),
        mean_residue = ('residue', 'mean'),
        std_residue = ('residue', 'std'),
        mean_solver = ('solver_total', 'mean'),
        std_solver = ('solver_total', 'std')
    )

    combined_cols = [
        hierarchy_data_summaries[['experiment', 'n_vertices', 'mean_hierarchy', 'std_hierarchy']],
        solver_ours_summaries[['median_iterations', 'mean_iterations', 'std_iterations', 'mean_solver', 'std_solver', 'mean_residue', 'std_residue']],
    ]
    if direct:
        combined_cols += [
            solver_direct_data[['direct_factor', 'direct_solve', 'pardiso_factor', 'pardiso_solve']]
        ]
    if sig21:
        combined_cols += [
            hierarchy_sig21_data[['sig21_hierarchy']],
            solver_sig21_data[['sig21_iterations', 'sig21_solver', 'sig21_residue']],
        ]
    if sig06:
        combined_cols += [
            hierarchy_sig06_data[['sig06_hierarchy']],
            solver_sig06_data[['sig06_iterations', 'sig06_solver', 'sig06_residue']],
        ]
    if amg:
        combined_cols += [
            amg_rs_data[['rs_hierarchy', 'rs_iterations', 'rs_solver']],
            amg_sa_data[['sa_hierarchy', 'sa_iterations', 'sa_solver']]
        ]
    combined_table = pd.concat(combined_cols, axis=1)
    combined_table = combined_table.sort_values('n_vertices')
    combined_table = combined_table.convert_dtypes()
    combined_table['experiment'] = combined_table['experiment'].replace('_', ' ', regex=True).str.title()
    combined_table['n_vertices'] = (combined_table['n_vertices'] / 1000).astype(int).astype(str) + 'k'
    combined_table['mean_hierarchy'] = combined_table['mean_hierarchy'] / 1000
    combined_table['std_hierarchy'] = combined_table['std_hierarchy'] / 1000
    combined_table['mean_solver'] = combined_table['mean_solver'] / 1000
    combined_table['std_solver'] = combined_table['std_solver'] / 1000
    combined_table['median_iterations'] = combined_table['median_iterations'].astype(int)
    combined_table['our_hierarchy'] = float_to_string(combined_table['mean_hierarchy'])
    combined_table['our_iterations'] = float_to_string(combined_table['mean_iterations'])
    combined_table['our_solve'] = float_to_string(combined_table['mean_solver'])
    combined_table['our_residue'] = float_to_string(combined_table['mean_residue'])
    if std:
        combined_table['our_hierarchy'] += '(' + float_to_string(combined_table['std_hierarchy']) + ')'
        combined_table['our_iterations'] += '(' + float_to_string(combined_table['std_iterations']) + ')'
        combined_table['our_solve'] += '(' + float_to_string(combined_table['std_solver']) + ')'
        combined_table['our_residue'] += '(' + float_to_string(combined_table['std_residue']) + ')'

    if direct:
        combined_table['direct_factor'] = combined_table['direct_factor'] / 1000
        combined_table['direct_solve'] = combined_table['direct_solve'] / 1000
        combined_table['pardiso_factor'] = combined_table['pardiso_factor'] / 1000
        combined_table['pardiso_solve'] = combined_table['pardiso_solve'] / 1000

    if sig21:
        combined_table['sig21_hierarchy'] = combined_table['sig21_hierarchy'] / 1000
        combined_table['sig21_solver'] = combined_table['sig21_solver'] / 1000

    if sig06:
        combined_table['sig06_hierarchy'] = combined_table['sig06_hierarchy'] / 1000
        combined_table['sig06_solver'] = combined_table['sig06_solver'] / 1000

    with pd.ExcelWriter(f'{out_dir}/{label}_{tau}_table.xlsx') as writer:
        combined_table.to_excel(writer)

    if (latex):
        latex_cols = []
        latex_headers = []
        if names_counts:
            latex_cols += ['experiment', 'n_vertices']
            latex_headers += ["Model", "Vertices"]
        latex_cols += ['our_hierarchy', 'median_iterations', 'our_solve']
        latex_headers = ["Hier. (s)", '#Iter.', 'Solve (s)']
        if sig21:
            latex_cols += ['sig21_hierarchy', 'sig21_iterations', 'sig21_solver']
            latex_headers += ['SIG21 Hier. (s)', '#Iter.', 'Solve (s)']
        if sig06:
            latex_cols += ['sig06_hierarchy', 'sig06_iterations', 'sig06_solver']
            latex_headers += ['SIG06 Hier. (s)', '#Iter.', 'Solve (s)']
        if amg:
            latex_cols += ['rs_hierarchy', 'rs_iterations', 'rs_solver', 'sa_hierarchy', 'sa_iterations', 'sa_solver']
            latex_headers += ['RS Hier. (s)', '#Iter.', 'Solve (s)', 'SA Hier. (s)', '#Iter.', 'Solve (s)']
        if direct:
            latex_cols += ['direct_factor', 'direct_solve', 'pardiso_factor', 'pardiso_solve']
            latex_headers += ['Eig. Fact. (s)', 'Eig. Subst. (s)', 'Par. Fact. (s)', 'Par. Subst. (s)']
        comparison_latex = tabulate(combined_table[latex_cols], headers=latex_headers,
        tablefmt="latex_booktabs", showindex="never", floatfmt=".2f")

        with (Path(out_dir).parents[0] / f'latex/comparisons_{label}_{tau}.tex').open('w') as f:
            f.write(comparison_latex)