import scanpy as sc
import squidpy as sq
import pandas as pd
import time
import os

from .utils import build_adata


def get_data(dataset_func=None, top_genes=3000, **args):
    start_time = time.time()
    assert (dataset_func), '>>> ERROR: You must appoint a function!'

    adata, n_cluster, dataset_details = dataset_func(**args)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print('>>> INFO: dataset name: {}, size: ({}, {}), cluster: {}.({:.3f}s)'.format(dataset_details, adata.X.shape[0], adata.X.shape[1], n_cluster, time.time() - start_time))

    return adata, n_cluster


def get_dlpfc_data(id, path=None):
    section_list = ['151507', '151508', '151509', '151510', '151669', '151670', 
                    '151671', '151672', '151673', '151674', '151675', '151676']
    if (isinstance(id, int) and id in range(12)):
        section_id = section_list[id]
    elif (isinstance(id, str) and id in section_list):
        section_id = id

    if (not path):
        path = os.path.join('st_datasets', 'DLPFC', section_id)
    adata = sc.read_visium(path=path)
    adata.var_names_make_unique()
    Ann_df = pd.read_csv(os.path.join(path, 'ground_truth.txt'), 
                        sep='\t', header=None, index_col=0)
    Ann_df.columns = ['Ground Truth']
    adata.obs['cluster'] = Ann_df.loc[adata.obs_names, 'Ground Truth']
    cluster_num = len(set(adata.obs['cluster'])) - 1
    return adata, cluster_num, 'dlpfc, slice id: {}'.format(section_id)


def get_human_breast_cancer_data(path=None):
    if (not path):
        path = 'st_datasets/human_breast_cancer/'
    adata = sc.read_visium(path=path)
    adata.var_names_make_unique()
    Ann_df = pd.read_csv(os.path.join(path, 'metadata.csv'), 
                        sep=',', header=None, index_col=0)
    Ann_df.columns = ['Ground Truth']
    adata.obs['cluster'] = Ann_df.loc[adata.obs_names, 'Ground Truth']
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, 'human_breast_cancer'


def get_mouse_brain_cerebellum_data(path=None):
    if (not path):
        path = 'st_datasets/mouse_brain_cerebellum/'
    adata = build_adata(spot_scale=1., 
                        root_path=path, 
                        name='mouse_brain_cerebellum')
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, 'mouse_brain_cerebellum'


def get_mouse_brain_sagittal_data(path=None, section=None, cut=None):
    assert (path and section and cut) or (not path and section and cut), \
        ">>> ERROR: Ambiguous appointment!"

    if (not path):
        path = os.path.join('st_datasets', 'mouse_brain_sagittal', 'section_'+str(section), cut)
    adata = build_adata(spot_scale=3., 
                        root_path=path, 
                        name='mouse_brain_sagittal')
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, 'mouse_brain_cerebellum'


def get_mouse_kidney_coronal_data(path=None):
    if (not path):
        path = 'st_datasets/mouse_kidney_coronal/'
    adata = build_adata(spot_scale=1., 
                        root_path=path, 
                        name='mouse_kidney_coronal')
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, 'mouse_kidney_coronal'


def get_mouse_olfactory_bulb_ST_data(id, path=None):
    if (not path):
        path = 'st_datasets/mouse_olfactory_bulb_ST/'
    adata = build_adata(spot_scale=10., 
                        root_path=os.path.join(path, 'slice_{}'.format(id)), 
                        name='mouse_olfactory_bulb_ST_slice_{}'.format(id))
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, 'mouse_olfactory_bulb_ST_slice_{}'.format(id)


def get_mouse_olfactory_bulb_Stereo_seq_data(path=None):
    if (not path):
        path = 'st_datasets/mouse_olfactory_bulb_Stereo-seq/'
    adata = build_adata(spot_scale=2., 
                        root_path=path, 
                        name='mouse_olfactory_bulb_Stereo-seq')
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, 'mouse_olfactory_bulb_Stereo-seq'


def get_mouse_somatosensory_cortex_data(path=None):
    if (not path):
        path = 'st_datasets/mouse_somatosensory_cortex/'
    adata = build_adata(spot_scale=2., 
                        root_path=path, 
                        name='mouse_somatosensory_cortex')
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, 'mouse_somatosensory_cortex'


def get_visium_hne_data(path=None):
    adata = sq.datasets.visium_hne_adata(path=path)
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, '10x_visium_hne'


def get_visium_fluo_data(path=None):
    adata = sq.datasets.visium_fluo_adata(path=path)
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster']))
    return adata, cluster_num, '10x_visium_fluo'
