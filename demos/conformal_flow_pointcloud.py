from pathlib import Path
# Base libraries
from plyfile import PlyData
import numpy as np
from robust_laplacian import point_cloud_laplacian
from scipy import sparse
# Viewer
import polyscope as ps
import polyscope.imgui as psim
# Multigrid Solver
from gravomg import MultigridSolver
from gravomg.util import neighbors_from_stiffness, normalize_bounding_box, normalize_axes

# Read point cloud
plydata = PlyData.read('../data/pointcloud/oil_pump.ply')
V = np.vstack((
    plydata['vertex']['x'],
    plydata['vertex']['y'],
    plydata['vertex']['z']
)).T

V = normalize_axes(V)

# Normalize bounding box and center
V = normalize_bounding_box(V)

# Compute operators
S, M = point_cloud_laplacian(V)
Minv = sparse.diags(1 / M.diagonal())

neigh = neighbors_from_stiffness(S)

# Camera parameters for screenshot script
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
ui_screenshot_dir = "/home/ruben/Dev/MultigridBilaplacianSolver/out/bottle/screenshots"

def conformal_flow(tau, n_iter, screenshot, screenshot_dir):
    global V, camera_pos
    path = Path(screenshot_dir)
    if screenshot:
        ps.screenshot(str(path / "0.jpg"))
    for i in range(n_iter):
        _, M = robust_laplacian.point_cloud_laplacian(V)
        lhs = M + tau * S
        rhs = M @ V
        V = solver.solve(lhs, rhs)
        V = normalize_bounding_box(V)
        pointcloud.update_point_positions(V)
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

pointcloud = ps.register_point_cloud('Input Point Cloud', V, enabled=True)
pointcloud.add_color_quantity('original positions', V, enabled=True)

ps.show()