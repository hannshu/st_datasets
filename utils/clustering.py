import scanpy as sc
import numpy as np
from sklearn.cluster import KMeans, Birch, AgglomerativeClustering
from sklearn.mixture import GaussianMixture


# mclust-EEE clustering method
# source: https://github.com/QIFEIDKN/STAGATE/blob/main/STAGATE/utils.py
def mclust_R(adata, num_cluster, modelNames='EEE', used_obsm='STAGATE', random_seed=0):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm[used_obsm]), num_cluster, modelNames)
    mclust_res = np.array(res[-2])

    return mclust_res.astype(int).astype(str)


def mclust(adata: sc.AnnData, n_cluster: int, use_rep: str = 'embedding') -> np.ndarray:
    return mclust_R(adata, used_obsm=use_rep, num_cluster=n_cluster)


def k_means(adata: sc.AnnData, n_cluster: int, use_rep: str = 'embedding') -> np.ndarray:
    return KMeans(n_clusters=n_cluster, random_state=0).fit_predict(adata.obsm[use_rep]).astype(str)


def birch(adata: sc.AnnData, n_cluster: int, use_rep: str = 'embedding') -> np.ndarray:
    return Birch(n_clusters=n_cluster).fit_predict(adata.obsm[use_rep]).astype(str)


def gmm(adata: sc.AnnData, n_cluster: int, use_rep: str = 'embedding') -> np.ndarray:
    return GaussianMixture(n_components=n_cluster, random_state=0).fit_predict(adata.obsm[use_rep]).astype(str)


def ahc(adata: sc.AnnData, n_cluster: int, use_rep: str = 'embedding') -> np.ndarray:
    return AgglomerativeClustering(n_clusters=n_cluster).fit_predict(adata.obsm[use_rep]).astype(str)
