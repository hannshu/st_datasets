import pandas as pd
import scanpy as sc
import numpy as np
import networkx as nx
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
    print(f">>> INFO: Generate {edge_list.shape[0]} edges, {edge_list.shape[0] / adata.shape[0]:.3f} edges per spot.({time.time() - start:.3f}s)")

    return edge_list


def convert_edge_to_adj(edge_list, spot_num=None, dense=True):
    g = build_nx_graph(edge_list, spot_num)
    
    if (dense):
        return nx.adjacency_matrix(g).todense()
    else:
        return nx.adjacency_matrix(g)


def convert_adj_to_edge(adj: np.array):
    return np.nonzero(adj)


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


def k_hop_neighbor(edge_list, adj, k, spot_num=None):
    if (edge_list):
        adj = convert_edge_to_adj(edge_list, spot_num)

    return adj ** k


def build_nx_graph(edge_list, spot_num=None):
    edge_list = [(edge_list[0][i], edge_list[1][i]) for i in range(len(edge_list[0]))]

    if (spot_num):
        g = nx.empty_graph(spot_num)
        g.add_edges_from(edge_list)
    else:
        g = nx.Graph(edge_list)

    return g