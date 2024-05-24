from typing import Optional
import numpy as np


def get_deconv_metric(label: np.ndarray, pred: np.ndarray, label_name: Optional[list] = None) -> None:

    from scipy.spatial.distance import jensenshannon
    from sklearn.metrics import mean_squared_error
    from skimage.metrics import structural_similarity

    def correlation_coefficient(T1, T2):
        numerator = np.mean((T1 - T1.mean()) * (T2 - T2.mean()))
        denominator = T1.std() * T2.std()
        if denominator == 0:
            return 0
        else:
            result = numerator / denominator
            return result
        
    def norm(x):
        if (1 < (x.max() - x.min())):
            return (x - np.min(x)) / (np.max(x) - np.min(x))
        return x
        
    label = norm(label)
    pred = norm(pred)

    if (None == label_name):
        label_name = list(range(label.shape[1]))
    label_wise_js = {
        name: f"{elem:.3f}" 
        for name, elem in zip(label_name, jensenshannon(label, pred))
    }

    print(f'>>> PCC(↑): {correlation_coefficient(label, pred):.3f}')
    print(f'>>> SSIM(↑): {structural_similarity(label, pred, data_range=1):.3f}')
    print(f'>>> RMSE(↓): {mean_squared_error(label, pred):.3f}')
    print(f'>>> JS(↓): {label_wise_js}')
    print(f'>>> JS(mean(↓)): {jensenshannon(label, pred).mean():.3f}')
