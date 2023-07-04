import igl
from plyfile import PlyData, PlyElement
import numpy as np

def read_mesh(file):
    if file.match('*.obj'):
        V, _, _, F, _, _ = igl.read_obj(str(file))
    elif file.match('*.off'):
        V, F, _ = igl.read_off(str(file))
    elif file.match('*.ply'):
        with file.open(mode='rb') as f:
            plydata = PlyData.read(f)
        V = np.stack([plydata['vertex'][el] for el in ['x', 'y', 'z']], axis=-1)
        F = np.stack(plydata['face']['vertex_indices'], axis=0).astype(np.int32)
    return V, F

def read_pointcloud(file):
    if not file.match('*.ply'): print('This function expects pointclouds in .ply format')
    assert file.match('*.ply')

    plydata = PlyData.read(str(file))  
    V = np.vstack((
        plydata['vertex']['x'],
        plydata['vertex']['y'],
        plydata['vertex']['z']
    )).T
    return V

def write_ply(pos, path):
    pos = pos.astype(np.float32)
    pos.dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
    pos = PlyElement.describe(pos.squeeze(), 'vertex')
    PlyData([pos]).write(path)
