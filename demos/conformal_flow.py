from pathlib import Path
import time
# Base libraries
import igl
import numpy as np
from scipy import sparse
# Viewer
import polyscope as ps
import polyscope.imgui as psim
# Multigrid Solver
from gravomg import MultigridSolver
from gravomg.util import neighbors_from_stiffness, normalize_area

from robust_laplacian import mesh_laplacian
from util import read_mesh

# Read mesh
V, F = read_mesh(Path('../data/nonmanifold/indonesian_statue.obj'))
 
# Normalize area and center around mean
V = normalize_area(V, F)

# Compute operators
M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_BARYCENTRIC)
S = -igl.cotmatrix(V, F)

Minv = sparse.diags(1 / M.diagonal())
neigh = neighbors_from_stiffness(S)

# Animate camera in screenshot script
camera_pos = np.array([-1, 1, -1])
angle = 0.002 * 2 * np.pi
cos, sin = np.cos(angle), np.sin(angle)
R = np.array([
    [cos, 0, -sin],
    [0, 1, 0],
    [sin, 0, cos]
])

# Create reusable solver
solver = MultigridSolver(V, neigh, M, tolerance=1e-4)

# Start UI
ui_tau = 0.01
ui_n_iter = 1
ui_screenshot = False
ui_screenshot_dir = "../out/screenshots"

def conformal_flow(tau, n_iter, screenshot, screenshot_dir):
    global V, camera_pos
    path = Path(screenshot_dir)
    if screenshot:
        ps.screenshot(str(path / "0.jpg"))
    for i in range(n_iter):
        M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_BARYCENTRIC)
        lhs = M + tau * S
        rhs = M @ V
        V = solver.solve(lhs, rhs)
        V = normalize_area(V, F)
        mesh.update_vertex_positions(V)
        if screenshot:
            camera_pos = R @ camera_pos
            ps.look_at(camera_pos, (0., 0., 0.))
            ps.screenshot(str(path / (str(i + 1) + ".jpg")))


# GUI
def conformal_flow_panel():
    global ui_tau, ui_n_iter, ui_screenshot, ui_screenshot_dir

    psim.TextUnformatted("Conformal Flow")
    changed, ui_tau = psim.InputFloat("tau", ui_tau, format='%.6f')
    changed, ui_n_iter = psim.InputInt("number of iterations", ui_n_iter)
    changed, ui_screenshot = psim.Checkbox("screenshot", ui_screenshot)
    changed, ui_screenshot_dir = psim.InputText("screenshot directory", ui_screenshot_dir)

    if(psim.Button("Compute")):
        conformal_flow(ui_tau, ui_n_iter, ui_screenshot, ui_screenshot_dir)

ps.init()
ps.set_ground_plane_mode('none')
ps.set_user_callback(conformal_flow_panel)

ps.look_at(camera_pos, (0., 0., 0.))

mesh = ps.register_surface_mesh('Input Mesh', V, F, smooth_shade=True, enabled=True)
mesh.add_color_quantity('normals', igl.per_vertex_normals(V, F), enabled=True)

ps.show()