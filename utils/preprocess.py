import pandas as pd
import scanpy as sc
import numpy as np
import igraph as ig
from sklearn.neighbors import NearestNeighbors
import time


def build_graph(adata: sc.AnnData, radius=None, knears=None, distance_metrics='l2', use_repo='spatial', matrix=None):
    start = time.time()

    if (isinstance(matrix, np.ndarray)):
        coor = pd.DataFrame(matrix)
    elif ('X' == use_repo):
        coor = pd.DataFrame(adata.X.todense())
    else:
        coor = pd.DataFrame(adata.obsm[use_repo])
        coor.index = adata.obs.index
        coor.columns = ['row', 'col']

    if (radius):
        nbrs = NearestNeighbors(radius=radius, metric=distance_metrics).fit(coor)
        _, indices = nbrs.radius_neighbors(coor, return_distance=True)
    else:
        nbrs = NearestNeighbors(n_neighbors=knears+1, metric=distance_metrics).fit(coor)
        _, indices = nbrs.kneighbors(coor)

    edge_list = np.array([[i, j] for i, sublist in enumerate(indices) for j in sublist])
    print(f">>> INFO: Generate {edge_list.shape[0]} edges, {(edge_list.shape[0] / adata.shape[0]) - 1:.3f} edges per spot.({time.time() - start:.3f}s)")

    return edge_list


def convert_edge_to_adj(edge_list, spot_num=None, dense=True):
    g = ig.Graph(n=spot_num, edges=edge_list.T)
    
    if (dense):
        return np.asarray(g.get_adjacency_sparse().todense())
    else:
        return g.get_adjacency_sparse()


def convert_adj_to_edge(adj: np.array):
    return np.array(np.nonzero(adj))


def get_not_adjacency_pair(adj: np.array):
    return np.array(np.nonzero(adj == 0))


def add_self_loop(edge_list, spot_num=None):
    if (spot_num):
        self_loop = np.arange(0, spot_num)
    else:
        self_loop = np.arange(0, np.max(edge_list))
    no_loop_list = delete_self_loop(edge_list)
    result_list = np.array([
        np.concatenate(no_loop_list[0], self_loop),
        np.concatenate(no_loop_list[1], self_loop)
    ])

    return result_list


def delete_self_loop(edge_list):
    non_duplicate_element = edge_list[0, :] != edge_list[1, :]
    filtered_array = np.array([
        edge_list[0][non_duplicate_element], 
        edge_list[1][non_duplicate_element]
    ])

    return filtered_array


def concat_adjacency_matrix(adata_list, edge_list, return_type=None):
    edges = edge_list[0]
    spot_count = adata_list[0].X.shape[0]

    for i in range(1, len(adata_list)):
        edge_list[i] = spot_count + edge_list[i]
        edges = np.vstack((edges, edge_list[i]))
        spot_count += adata_list[i].X.shape[0]

    if ('adj' == return_type):
        return convert_edge_to_adj(edges.T)
    else:
        return edges.T
    

def conv_to_one_hot(X, n):
    X = np.array(X).astype(int)
    refer = np.eye(n)
    X_one_hot = refer[X]
    return X_one_hot


def k_hop_adj(edge_list, adj, k, spot_num=None):
    if (edge_list):
        adj = convert_edge_to_adj(edge_list, spot_num)

    return adj ** k
