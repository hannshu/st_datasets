import numpy as np
from sklearn.preprocessing import LabelEncoder
from rpy2.rinterface_lib.embedded import RRuntimeError


# mclust-EEE clustering method
# source: https://github.com/QIFEIDKN/STAGATE/blob/main/STAGATE/utils.py
def mclust_R(adata, num_cluster, modelNames='EEE', used_obsm='STAGATE', random_seed=2020):
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

    adata.obs['mclust'] = mclust_res
    adata.obs['mclust'] = adata.obs['mclust'].astype('int')
    adata.obs['mclust'] = adata.obs['mclust'].astype('category')
    return adata


def evaluate_embedding(adata, n_cluster, cluster_method=['mclust'], cluster_score_method='ARI'):
    # evaluate cluster score
    obs_df = adata.obs.dropna()
    true_label = LabelEncoder().fit_transform(obs_df['cluster'])
    score = {}

    try:
        if ('mclust' in cluster_method):
            adata = mclust_R(adata, used_obsm='embedding', num_cluster=n_cluster)
            obs_df = adata.obs.dropna()
            pred_label = LabelEncoder().fit_transform(obs_df['mclust']).astype(str)
            score['mclust'] = cal_cluster_score(true_label, pred_label, cluster_score_method)
            print(f">>> INFO: Finish mclust clustering, clustering score {score['mclust']:.3f}")
    except (TypeError or RRuntimeError):
        print('>>> WARNING: mclust report TypeError')
        score['mclust'] = -1

    if ('kmeans' in cluster_method):
        from sklearn.cluster import KMeans
        adata.obs['kmeans'] = KMeans(n_clusters=n_cluster, random_state=0).fit_predict(adata.obsm['embedding']).astype(str)
        obs_df = adata.obs.dropna()
        pred_label = LabelEncoder().fit_transform(obs_df['kmeans'])
        score['kmeans'] = cal_cluster_score(true_label, pred_label, cluster_score_method)
        print(f">>> INFO: Finish KMeans clustering, clustering score {score['kmeans']:.3f}")

    if ('birch' in cluster_method):
        from sklearn.cluster import Birch
        adata.obs['birch'] = Birch(n_clusters=n_cluster).fit_predict(adata.obsm['embedding']).astype(str)
        obs_df = adata.obs.dropna()
        pred_label = LabelEncoder().fit_transform(obs_df['birch'])
        score['birch'] = cal_cluster_score(true_label, pred_label, cluster_score_method)
        print(f">>> INFO: Finish birch clustering, clustering score {score['birch']:.3f}")

    if ('gmm' in cluster_method):
        from sklearn.mixture import GaussianMixture
        adata.obs['gmm'] = GaussianMixture(n_components=adata.obsm['embedding'].shape[1], random_state=0).fit_predict(adata.obsm['embedding']).astype(str)
        obs_df = adata.obs.dropna()
        pred_label = LabelEncoder().fit_transform(obs_df['gmm'])
        score['gmm'] = cal_cluster_score(true_label, pred_label, cluster_score_method)
        print(f">>> INFO: Finish GMM clustering, clustering score {score['gmm']:.3f}")

    if ('ahc' in cluster_method):
        from sklearn.cluster import AgglomerativeClustering   
        adata.obs['ahc'] = AgglomerativeClustering(n_clusters=n_cluster).fit_predict(adata.obsm['embedding']).astype(str)
        obs_df = adata.obs.dropna()
        pred_label = LabelEncoder().fit_transform(obs_df['ahc'])
        score['ahc'] = cal_cluster_score(true_label, pred_label, cluster_score_method)
        print(f">>> INFO: Finish AHC clustering, clustering score {score['ahc']:.3f}")

    return adata, score


def cal_cluster_score(true_label, pred_label, method='ARI'):
    if ('ARI' == method):
        from sklearn.metrics import adjusted_rand_score
        return adjusted_rand_score(true_label, pred_label)
    elif ('NMI' == method):
        from sklearn.metrics import normalized_mutual_info_score
        return normalized_mutual_info_score(true_label, pred_label)
    elif ('AMI' == method):
        from sklearn.metrics import adjusted_mutual_info_score 
        return adjusted_mutual_info_score(true_label, pred_label)
    elif ('v_measure_score' == method):
        from sklearn.metrics import v_measure_score
        return v_measure_score(true_label, pred_label)
    elif ('silhouette_score' == method):
        from sklearn.metrics import silhouette_score
        # if counting silhouette_score, you have to use adata[obs_df.rows] as true_label
        return silhouette_score(true_label, pred_label)
    else:
        print(f'>>> ERROR: Method {method} is not supported!')
