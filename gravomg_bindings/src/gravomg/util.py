import numpy as np
from scipy.spatial import KDTree

def neighbors_from_stiffness(S):
    S_coo = S.tocoo()
    node_j, node_i = S_coo.row, S_coo.col
    
    return homogenize_edges(node_i, node_j)

def neighbors_from_faces(F):
    node_i = np.concatenate([F[:, 0], F[:, 0], F[:, 1], F[:, 1], F[:, 2], F[:, 2]])
    node_j = np.concatenate([F[:, 1], F[:, 2], F[:, 0], F[:, 2], F[:, 0], F[:, 1]])
    return homogenize_edges(*coalesce_edges(node_i, node_j))

def knn(V, k):
    tree = KDTree(V)
    return tree.query(V, k + 1)[1][:, 1:]

def knn_undirected(V, k):
    node_i = np.arange(V.shape[0])[:, None].repeat(k, axis=1).flatten()
    node_j = knn(V, k).flatten()

    # Make undirected
    node_i = np.concatenate([node_i, node_j])
    node_j = np.concatenate([node_j, node_i])

    return homogenize_edges(*coalesce_edges(node_i, node_j))

def coalesce_edges(node_i, node_j):
    sort_idx = np.argsort(node_i) 
    node_i = node_i[sort_idx]
    node_j = node_j[sort_idx]
    node_i, node_j = np.unique(np.stack([node_i, node_j], axis=0), axis=1)
    return node_i, node_j

def homogenize_edges(node_i, node_j):
    nodes, degree = np.unique(node_i, return_counts=True)
    k_new = degree.max()

    # Homogenize indices
    padded_idx = np.arange(node_j.shape[0]) - np.pad(np.cumsum(degree), (1, 0))[node_i] + node_i * k_new
    neigh = -1 * np.ones(nodes.shape[0] * k_new, dtype=np.int32)
    neigh[padded_idx] = node_j
    return neigh.reshape(-1, k_new)

def face_area(pos, F):
    v1 = pos[F[:, 0]]
    v2 = pos[F[:, 1]]
    v3 = pos[F[:, 2]]
    return np.linalg.norm(np.cross(v2 - v1, v3 - v1), axis=1) / 2

def normalize_area(pos, F):
    pos = pos / np.sqrt(face_area(pos, F).sum());
    pos = pos - np.mean(pos, axis=0, keepdims=True)
    return pos

def normalize_bounding_box(pos):
    pos = pos - pos.mean(axis=0, keepdims=True)
    scale = 0.5 / np.abs(pos).max()
    return pos * scale

def normalize_axes(pos):
    _, u, v = np.linalg.svd(pos)
    v = (v.T[np.argsort(u)]).T
    return (v @ pos.T).T

def normalize_axes(pos):
    ax_permutation = np.argsort(np.std(pos, axis=0))
    return pos[:, ax_permutation]