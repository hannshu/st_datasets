import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import os
import time


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