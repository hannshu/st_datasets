import scanpy as sc
import numpy as np
from math import floor
from typing import Optional, Union, List


def get_micro_env_img(
    adata: sc.AnnData, 
    offset: Union[int, float],
    library_id: Optional[str] = None, 
    img_key: str = "hires",
    scale_factor: Optional[float] = None,
    spatial_key: str = 'spatial'
) -> List[np.ndarray]:

    # get image, coordinates and scale factor
    # img.shape: (y axis, x axis, RGB)
    # coordinates.shape: (cells, (x, y))
    if (library_id):
        _, spatial_data = sc.pl._tools.scatterplots._check_spatial_data(adata.uns, library_id)
    else:
        _, spatial_data = sc.pl._tools.scatterplots._check_spatial_data(adata.uns, sc._utils._empty)
    img, img_key = sc.pl._tools.scatterplots._check_img(spatial_data, img=None, img_key=img_key)
    scale_factor = sc.pl._tools.scatterplots._check_scale_factor(
        spatial_data, img_key=img_key, scale_factor=scale_factor
    )
    
    # add padding for the image to make sure 
    # 1. each block has the same shape
    # 2. cell always located at the center of the block
    img = np.pad(img, ((offset, offset), (offset, offset), (0, 0)), 'constant', constant_values=(1.,))   
    coords = adata.obsm[spatial_key] * scale_factor # fit coords to the image

    return [
        img[
            floor(coor[1]): floor(coor[1]+2*offset),    # for y axis
            floor(coor[0]): floor(coor[0]+2*offset),    # for x axis
            :
        ] for coor in coords
    ]   # return a matrix with each cell's micro env image, shaped (2*offset, 2*offset)
