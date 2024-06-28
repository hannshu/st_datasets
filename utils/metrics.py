import numpy as np
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    adjusted_mutual_info_score,
    v_measure_score,
    silhouette_score,
    davies_bouldin_score
)
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import mean_squared_error
from skimage.metrics import structural_similarity


# clustering metrics
def cal_ARI(label: np.ndarray, pred: np.ndarray) -> float:
    return adjusted_rand_score(label['nan' != label], pred['nan' != label])


def cal_NMI(label: np.ndarray, pred: np.ndarray) -> float:
    return normalized_mutual_info_score(label['nan' != label], pred['nan' != label])


def cal_AMI(label: np.ndarray, pred: np.ndarray) -> float:
    return adjusted_mutual_info_score(label['nan' != label], pred['nan' != label])


def cal_v_measure_score(label: np.ndarray, pred: np.ndarray) -> float:
    return v_measure_score(label['nan' != label], pred['nan' != label])


def cal_SC(embedding: np.ndarray, pred: np.ndarray) -> float:
    return silhouette_score(embedding, pred)


def cal_DB(embedding: np.ndarray, pred: np.ndarray) -> float:
    return davies_bouldin_score(embedding, pred)


def norm(x: np.ndarray) -> np.ndarray:
    if (1 < (x.max() - x.min())):
        return (x - np.min(x)) / (np.max(x) - np.min(x))
    return x


# deconvolution metrics
def cal_PCC(label: np.ndarray, pred: np.ndarray) -> float:
    return pearsonr(norm(label).flatten(), norm(pred).flatten())[0]


def cal_SPCC(label: np.ndarray, pred: np.ndarray) -> float:
    return spearmanr(norm(label).flatten(), norm(pred).flatten())[0]


def cal_SSIM(label: np.ndarray, pred: np.ndarray) -> float:
    return structural_similarity(norm(label), norm(pred), data_range=1)


def cal_RMSE(label: np.ndarray, pred: np.ndarray) -> float:
    return np.sqrt(mean_squared_error(norm(label), norm(pred)))


def cal_JSD(label: np.ndarray, pred: np.ndarray) -> float:
    return jensenshannon(norm(label), norm(pred)).mean()
