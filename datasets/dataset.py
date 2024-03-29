import scanpy as sc
import numpy as np
import time
import os

from .utils import check_file_location


home_dir = os.environ['HOME']
dataset_url = 'https://huggingface.co/datasets/han-shu/st_datasets/resolve/main'


def datasets():
    print(">>> INFO: find available datasets:")
    for key, value in globals().items(): 
        if (callable(value) and key.startswith('get_') and 'get_data' != key): 
            print(value)


def get_data(dataset_func=None, top_genes=3000, preprocess=True, **args):
    start_time = time.time()
    assert (dataset_func), '>>> ERROR: You must appoint a function!'

    if (dataset_func in globals().values()):
        adata, n_cluster, dataset_details = dataset_func(**args)
    else:
        adata = dataset_func(**args)
        n_cluster = None
        dataset_details = str(dataset_func)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=top_genes)
    if (preprocess):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    if ('cluster' in adata.obs):
        adata.obs['cluster'] = adata.obs['cluster'].astype(str)
    print('>>> INFO: dataset name: {}, size: ({}, {}), cluster: {}.({:.3f}s)'.format(
        dataset_details, adata.X.shape[0], adata.X.shape[1], n_cluster, time.time() - start_time
    ))

    return adata, n_cluster


def get_dlpfc_data(id, path=None):
    section_list = ['151507', '151508', '151509', '151510', '151669', '151670', 
                    '151671', '151672', '151673', '151674', '151675', '151676']

    if (isinstance(id, int) and id in range(12)):
        section_id = section_list[id]
    elif (isinstance(id, str) and id in section_list):
        section_id = id

    if (not path):
        path = os.path.join(home_dir, 'data', 'DLPFC', f'DLPFC_{section_id}.h5ad')
    url = f'{dataset_url}/DLPFC/DLPFC_{section_id}.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster'])) - 1

    return adata, cluster_num, f'dorsolateral prefrontal cortex (DLPFC), slice: {section_id}'


def get_human_breast_cancer_data(path=None):
    if (not path):
        path = os.path.join(home_dir, 'data', 'human_breast_cancer.h5ad')
    url = f'{dataset_url}/human_breast_cancer.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, 'human breast cancer'


def get_mouse_brain_ffpe_data(path=None):
    if (not path):
        path = os.path.join(home_dir, 'data', 'FFPE_Mouse_Brain.h5ad')
    url = f'{dataset_url}/FFPE_Mouse_Brain.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, 'mouse brain coronal FFPE data'


def get_mouse_brain_sagittal_data(path=None, section: int=None, pos: str=None):
    assert (section in [1, 2]) and (pos in ['Ant', 'Pos']), ">>> ERROR: Ambiguous appointment!"
    
    file_list = [
        {
            'Ant': ['STDS0000018', 'Mouse_Brain_Serial_Section_1_Sagittal-Anterior_10xvisium_processed.h5ad'],
            'Pos': ['STDS0000021', 'Mouse_Brain_Serial_Section_1_Sagittal-Posterior_10xvisium_processed.h5ad'],
        }, 
        {
            'Ant': ['STDS0000019', 'Mouse_Brain_Serial_Section_2_Sagittal-Anterior_10xvisium_processed.h5ad'],
            'Pos': ['STDS0000020', 'Mouse_Brain_Serial_Section_2_Sagittal-Posterior_10xvisium_processed.h5ad'],
        }
    ] 

    if (not path):
        path = os.path.join(home_dir, 'data', 'MBS', file_list[section-1][pos][1])
    url = f'https://ftp.cngb.org/pub/SciRAID/stomics/{file_list[section-1][pos][0]}/stomics/{file_list[section-1][pos][1]}'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    adata.obs['cluster'] = adata.obs['clusters']
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, f'mouse brain sagittal, section: {section}, position: {pos}'


def get_mouse_brain_cerebellum_data(path=None):
    if (not path):
        path = os.path.join(home_dir, 'data', 'mouse_brain_cerebellum.h5ad')
    url = f'{dataset_url}/mouse_brain_cerebellum.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, 'mouse brain cerebellum'


def get_mouse_kidney_coronal_data(path=None):
    if (not path):
        path = os.path.join(home_dir, 'data', 'mouse_kidney_coronal.h5ad')
    url = f'{dataset_url}/mouse_kidney_coronal.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, 'mouse kidney coronal'


def get_mouse_olfactory_bulb_data(path=None, tech=None, **args):
    assert (tech in ['visium', 'Stereo-seq', 'ST', 'Slide_seqV2']), ">>> ERROR: Ambiguous appointment!"
    
    if ('ST' == tech):
        assert ('id' in args), ">>> ERROR: Ambiguous appointment!"
        section_list = [
            'Rep1_MOB_ST_deal.h5ad', 'Rep2_MOB_ST_deal.h5ad', 'Rep3_MOB_ST_deal.h5ad', 'Rep4_MOB_ST_deal.h5ad',
            'Rep5_MOB_ST_deal.h5ad', 'Rep6_MOB_ST_deal.h5ad', 'Rep7_MOB_ST_deal.h5ad', 'Rep8_MOB_ST_deal.h5ad',
            'Rep9_MOB_ST_deal.h5ad', 'Rep10_MOB_ST_deal.h5ad', 'Rep11_MOB_ST_deal.h5ad', 'Rep12_MOB_ST_deal.h5ad'
        ]
        if (not path):
            path = os.path.join(home_dir, 'data', 'mouse_olfactory_bulb', section_list[args['id']])
        url = f'{dataset_url}/mouse_olfactory_bulb/ST/{section_list[args["id"]]}'
    elif ('Stereo-seq' == tech):
        if (not path):
            path = os.path.join(home_dir, 'data', 'mouse_olfactory_bulb', 'mob_stereo_seq.h5ad')
        url = f'{dataset_url}/mouse_olfactory_bulb/mob_stereo_seq.h5ad'
    elif ('Slide_seqV2' == tech):
        if (not path):
            path = os.path.join(home_dir, 'data', 'mouse_olfactory_bulb', 'Mob_200127_15.h5ad')
        url = f'{dataset_url}/mouse_olfactory_bulb/Mob_200127_15.h5ad'
    else:
        if (not path):
            path = os.path.join(home_dir, 'data', 'mouse_olfactory_bulb', 'GSM4656181_10x_Visium_deal.h5ad')
        url = f'{dataset_url}/mouse_olfactory_bulb/GSM4656181_10x_Visium_deal.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()

    if (not (('Stereo-seq' == tech) or ('Slide_seqV2' == tech))):
        adata.obs['cluster'] = adata.obs['clusters']
        cluster_num = len(set(adata.obs['cluster']))
    else:
        cluster_num = None

    return adata, cluster_num, f'mouse olfactory bulb obtained by {tech}'


def get_mouse_somatosensory_cortex_data(path=None):
    if (not path):
        path = os.path.join(home_dir, 'data', 'mouse_somatosensory_cortex.h5ad')
    url = f'{dataset_url}/mouse_somatosensory_cortex.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, 'mouse somatosensory cortex'


def get_zesta_data(path=None, id=None):
    if (not path):
        path = os.path.join(home_dir, 'data', 'spatial_sixtime_slice_stereoseq.h5ad')
    url = 'https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000057/stomics/spatial_sixtime_slice_stereoseq.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()

    if (id in ['10hpf', '12hpf', '18hpf', '24hpf', '3.3hpf', '5.25hpf']):
        adata = adata[adata.obs['time'] == id, :]

    adata.obsm['spatial'] = np.array([np.array(adata.obs['spatial_x']), np.array(adata.obs['spatial_y'])]).T
    adata.obs['cluster'] = adata.obs['bin_annotation']
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, 'zebrafish embryogenesis spatiotemporal transcriptomic atlas (ZESTA)'


def get_mosta_data(path=None, time='E9.5', id='E1S1'):
    if (not path):
        path = os.path.join(home_dir, 'data', 'MOSTA', f'{time}_{id}.MOSTA.h5ad')
    url = f'https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/{time}_{id}.MOSTA.h5ad'

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()
    adata.obs['cluster'] = adata.obs['annotation']
    cluster_num = len(set(adata.obs['cluster']))

    return adata, cluster_num, f'mouse organogenesis spatiotemporal transcriptomic atlas (MOSTA), section: {id}'


def get_outside_data(path=None, url=None):
    assert (path) or (url), '>>> ERROR: at lease provide a prth or a url!'

    if (None == path):
        # if obtain data from remote, download to ~/data/user_downloaded/data_name.h5ad
        path = os.path.join(home_dir, 'data', 'user_downloaded', url.split('/')[-1])

    adata = sc.read_h5ad(check_file_location(path=path, url=url))
    adata.var_names_make_unique()

    return adata, None, f'Outside dataset'
