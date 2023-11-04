import pandas as pd
import scanpy as sc
import numpy as np
import networkx as nx
from sklearn.neighbors import NearestNeighbors


def build_graph(adata: sc.AnnData, radius=None, knears=None):
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['row', 'col']

    if (radius):
        nbrs = NearestNeighbors(radius=radius).fit(coor)
        _, indices = nbrs.radius_neighbors(coor, return_distance=True)
    else:
        nbrs = NearestNeighbors(n_neighbors=knears+1).fit(coor)
        _, indices = nbrs.kneighbors(coor)

    edge_list = []
    for i in range(len(indices)):
        for j in range(len(indices[i])):
            edge_list.append([i, indices[i][j]])

    return edge_list


def convert_edge_to_adj(edge_list, spot_num=None, dense=True):
    edge_list = [(edge_list[0][i], edge_list[1][i]) for i in range(len(edge_list[0]))]

    if (spot_num):
        g = nx.empty_graph(spot_num)
        g.add_edges_from(edge_list)
    else:
        g = nx.Graph(edge_list)
    
    if (dense):
        return nx.adjacency_matrix(g).todense()
    else :
        return nx.adjacency_matrix(g)


def convert_adj_to_edge(adj: np.array):
    return np.nonzero(adj)

