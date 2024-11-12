import pandas as pd
import scanpy as sc
import numpy as np
import igraph as ig
from typing import Optional
from scipy.sparse import csr_matrix, issparse
from sklearn.neighbors import NearestNeighbors
from collections import Counter
import time


def get_feature(adata: sc.AnnData, hvg: bool = False) -> np.ndarray:

    def get_x(adata):
        return adata.X.toarray() if (issparse(adata.X)) else adata.X

    if (False == hvg):
        return get_x(adata)
    else:
        return get_x(adata[:, adata.var.highly_variable])


# reference: https://github.com/BayraktarLab/cell2location/issues/204#issuecomment-1272416837
def inverse_lognormalize(adata: sc.AnnData, use_count: str = 'total_counts') -> sc.AnnData:

    adata = adata.copy()

    normalized_data = np.exp(get_feature(adata)) - 1
    data = normalized_data / 10000 * adata.obs[use_count].to_numpy().reshape(-1, 1)
    adata.X = csr_matrix(data.astype(int))

    return adata


def build_graph(adata, radius=None, knears=None, distance_metrics='l2', use_repo='spatial'):
    start = time.time()

    if (isinstance(adata, np.ndarray)):
        coor = pd.DataFrame(adata)
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
    

def conv_to_one_hot(X):
    from sklearn.preprocessing import OneHotEncoder

    enc = OneHotEncoder(sparse=False)
    one_hot = enc.fit_transform(X.reshape(-1, 1)) 
    return one_hot


def k_hop_adj(edge_list, adj, k, spot_num=None):
    if (edge_list):
        adj = convert_edge_to_adj(edge_list, spot_num)

    return adj ** k


def get_grid(coors: np.ndarray, distance: float) -> np.ndarray:

    x_scale = [coors[:, 0].min(), coors[:, 0].max()]
    y_scale = [coors[:, 1].min(), coors[:, 1].max()]

    return np.meshgrid(
        np.linspace(x_scale[0], x_scale[1], int((x_scale[1] - x_scale[0]) // distance)), 
        np.linspace(y_scale[0], y_scale[1], int((y_scale[1] - y_scale[0]) // distance))
    )


def get_simulated_data(
    adata: sc.AnnData, 
    distance: float, 
    label_tag: Optional[str] = 'cluster', 
    loc_style: str = 'visium', 
    used_obsm: str = 'spatial', 
    aggr_method: str = 'sum',
) -> sc.AnnData:

    start_time = time.time()
    assert(used_obsm in adata.obsm), '>>> ERROR: No coordinations in adata!'

    xs, ys = get_grid(adata.obsm[used_obsm], distance)
    array_row = np.tile(np.arange(xs.shape[0]), (xs.shape[1], 1)).T.reshape(-1)
    array_col = np.tile(np.arange(xs.shape[1]), (xs.shape[0], 1)).reshape(-1)

    if ('visium' == loc_style):
        xs[list(range(1, xs.shape[0], 2)), :] += distance / 2
    coors = np.array([xs.reshape(-1), ys.reshape(-1)]).T

    nbrs = NearestNeighbors(radius=distance/2).fit(adata.obsm[used_obsm])
    neighbors = nbrs.radius_neighbors(coors, return_distance=False)
    entity_item = np.array([len(item) != 0 for item in neighbors])
    
    neighbors = neighbors[entity_item]
    coors = coors[entity_item]
    array_row = array_row[entity_item]
    array_col = array_col[entity_item]

    assert (np.mean([len(neighbors[i]) for i in range(len(neighbors))]) > 1.0), \
        '>>> ERROR: Invalid result, try to increase `distance` parameter.'

    x = get_feature(adata)
    aggr_func = np.sum if ('sum' == aggr_method) else np.mean  
    x = np.array([aggr_func(x[neighbors[i]], axis=0) for i in range(len(neighbors))])

    obs_df = pd.DataFrame(
        data=np.array([array_row, array_col, [len(neighbors[i]) for i in range(len(neighbors))]]).T
        if (not 'total_counts' in adata.obs) else np.array([
            array_row, array_col, 
            [len(neighbors[i]) for i in range(len(neighbors))],
            [int(aggr_func(adata.obs['total_counts'][neighbors[i]], axis=0)) for i in range(len(neighbors))]
        ]).T, 
        index=[
            f'{list(adata.obs.index[neighbors[i]].astype(str))}'
            for i in range(len(neighbors))
        ], 
        columns=['array_row', 'array_col', 'cell_count', 'total_counts'] 
        if ('total_counts' in adata.obs) else ['array_row', 'array_col', 'cell_count']
    )

    sim_adata = sc.AnnData(csr_matrix(x), obs=obs_df, var=adata.var)
    sim_adata.obsm[used_obsm] = coors

    if (label_tag):
        freq = [Counter(adata.obs[label_tag][neighbors[i]]) for i in range(len(neighbors))]

        sim_adata.obs['cluster'] = [item.most_common(1)[0][0] for item in freq]
        sim_adata.obsm['type_count'] = pd.DataFrame(data=freq, index=sim_adata.obs.index).fillna(0).astype(int)
        sim_adata.obsm['deconvolution_result'] = sim_adata.obsm['type_count'].div(sim_adata.obsm['type_count'].sum(axis=1), axis=0)

        print('>>> INFO: Ground truth cell type proportion matrix is saved at sim_adata.obsm[\'deconvolution_result\']')

    print(f'>>> INFO: Generate simulated ST data, simulated data shape {sim_adata.shape}, average {np.mean([len(neighbors[i]) for i in range(len(neighbors))]):.3f} cells in each spot. ({time.time() - start_time:.3f}s)')

    return sim_adata
