import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import numpy as np
from scipy.sparse import coo_matrix


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


def visualize_graph(adata: sc.AnnData, edges, save_path=None, color_label=None):
    without_ring = np.array([[edge[0], edge[1]] for edge in edges if (edge[0] != edge[1])]).T
    if ((0, ) != without_ring.shape):
        adata.obsp['visualize'] = coo_matrix(
            ([1]*without_ring.shape[1], (without_ring)), 
            shape=(adata.X.shape[0], adata.X.shape[0]), 
            dtype=np.int32
        )
        
        if (None == color_label):
            adata.obs['default'] = ['1'] * adata.X.shape[0]
            color_label = ['default']

        sq.pl.spatial_scatter(
            adata, connectivity_key="visualize", img=False, color=color_label, legend_loc=None,
            edges_width=0.25, title='', axis_label=['', ''], edges_color='black'
        )
        if (save_path):
            plt.savefig(save_path, dpi=1600)

        del adata.obsp['visualize']
    else:
        print(">>> WARNING: No edge to draw!")

