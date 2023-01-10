import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import os
import time
from sklearn.neighbors import NearestNeighbors
import seaborn as sns


def get_data(dataset_func=None, top_genes=3000, **args):
    start_time = time.time()
    assert (dataset_func), '>>> ERROR: You must appoint a function!'

    adata, n_cluster, dataset_details = dataset_func(**args)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print('>>> INFO: dataset name: {}, size: ({}, {}), cluster: {}.({:.3f}s)'.format(
        dataset_details, adata.X.shape[0], adata.X.shape[1], n_cluster, time.time() - start_time))

    return adata, n_cluster


def build_adata(spot_scale, root_path=None, name='dataset'):
    assert root_path != None, ">>> ERROR: root_path cannot be None!"

    x_path = os.path.join(root_path, 'x.gz')
    true_path = os.path.join(root_path, 'ground_truth.csv')
    spatial_path = os.path.join(root_path, 'spatial.csv')
    img_path = os.path.join(root_path, 'hist.png')

    df_x = pd.read_csv(x_path, sep='\t')
    x = np.array(np.array(df_x).T[1:], dtype='float')
    x = sparse.csr_matrix(x)

    df_true = pd.read_csv(true_path)
    true = np.array(df_true).T[1]
    true_name = np.array(df_true).T[0]
    true = np.array([str(k) for k in true])

    spa = pd.read_csv(spatial_path)
    spa = spa.drop(labels="var_names", axis=1)
    spa = np.array(spa)

    img=plt.imread(img_path)

    adata = sc.AnnData(x)
    adata.var_names = df_x['var_names']
    adata.obs_names = true_name
    adata.obsm['spatial'] = np.array(spa)
    adata.uns['spatial'] = {}
    adata.uns['spatial'][name] = {}
    adata.uns['spatial'][name]['images'] = {}
    adata.uns['spatial'][name]['images']['hires'] = img
    adata.uns['spatial'][name]['scalefactors'] = {}
    adata.uns['spatial'][name]['scalefactors']['spot_diameter_fullres'] = spot_scale
    adata.uns['spatial'][name]['scalefactors']['tissue_hires_scalef'] = 1.
    # adata.uns['spatial'][name]['scalefactors']['fiducial_diameter_fullres'] = 1.
    # adata.uns['spatial'][name]['scalefactors']['tissue_lowres_scalef'] = 1.
    adata.obs['cluster'] = true

    return adata


def build_graph(adata, radius=None, knears=None, format='adj', save_img=None, **args):
    assert (radius and not knears) or (not radius and knears), \
        ">>> ERROR: You can only choose one method."
    assert format in ['adj', 'pyg'], \
        ">>> ERROR: This return data type is not provided now."

    if ('highly_variable' in adata.var):
        adata = adata[:, adata.var['highly_variable']]

    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['row', 'col']
    if (radius):
        nbrs = NearestNeighbors(radius=radius).fit(coor)
        _, indices = nbrs.radius_neighbors(coor, return_distance=True)
    else:
        nbrs = NearestNeighbors(n_neighbors=knears+1).fit(coor)
        _, indices = nbrs.kneighbors(coor)

    adj = np.zeros((adata.X.shape[0], adata.X.shape[0]), dtype='float')
    for i in range(len(indices)):
        for j in range(len(indices[i])):
            adj[i][indices[i][j]] = 1.

    if ('pyg' == format):
        import torch
        from torch_geometric.data import Data

        data = Data(edge_index=torch.LongTensor(np.nonzero(adj)), 
                    x=torch.FloatTensor(adata.X.todense()))

        if (save_img):
            visualize_graph(adata, data.edge_index.T, save_img, **args)

        return data
    else:
        if (save_img):
            visualize_graph(adata, np.array(np.nonzero(adj)).T, save_img, **args)

        return adj


def visualize_graph(adata: sc.AnnData, edges, save_path, **args):
    import squidpy as sq
    import torch
    from torch_geometric import utils

    edge_list_without_ring = np.array([[edge[0], edge[1]] for edge in edges if (edge[0] != edge[1])]).T
    adata.obsp['visualize'] = utils.to_scipy_sparse_matrix(torch.LongTensor(edge_list_without_ring))
    sq.pl.spatial_scatter(adata, connectivity_key="visualize", img=False, **args)
    plt.savefig(save_path, dpi=1600)
    del adata.obsp['visualize']


def show_distrib_map(adata, title='data_distribution_map', format='hist', **args):
    assert format in ['hist', 'violin'], \
        ">>> ERROR: This return data type is not provided now."

    if ('highly_variable' in adata.var):
        adata = adata[:, adata.var['highly_variable']]

    val = adata.X.todense().reshape(adata.X.shape[0]*adata.X.shape[1]).tolist()[0]
    if ('violin' == format):
        df = pd.DataFrame()
        df['gen_exp'] = val
        df['name'] = [title]*len(val)
        sns.violinplot(data=df, **args) 
        plt.title(title)
    else:
        sns.histplot(val, **args)
        plt.title(title)
    plt.savefig(title+'.png', dpi=600)
