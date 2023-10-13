# CONTENT 

| **dataset**                        | **technology**      | **slices** | **spots/cells** | **genes**     | **source**                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|------------------------------------|---------------------|------------|-----------------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DLPFC dataset                      | 10x Genomics Visium | 12         | 3460 - 4789     | 33538         | [link](https://doi.org/10.1038/s41593-020-00787-0)                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| human breast cancer dataset        | 10x Genomics Visium | 1          | 3798            | 36601         | [link](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_1 )                                                                                                                                                                                                                                                                                                                                                                               |
| mouse olfactory bulb dataset       | ST                  | 12         | 231 - 282       | 15284 - 16675 | [link](https://www.science.org/doi/10.1126/science.aaf2403)                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| mouse kidney coronal dataset       | 10x Genomics Visium | 1          | 1438            | 31053         | [link](https://www.10xgenomics.com/resources/datasets/mouse-kidney-section-coronal-1-standard-1-1-0)                                                                                                                                                                                                                                                                                                                                                                                             |
| mouse brain sagittal dataset       | 10x Genomics Visium | 4          | 2696 - 3353     | 31053         | [link](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0), [link](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-1-0), [link](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard-1-1-0), [link](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-posterior-1-standard-1-1-0) |
| mouse somatosensory cortex dataset | osmFISH             | 1          | 5328            | 33            | [link](https://doi.org/10.1038/s41592-018-0175-z)                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| mouse olfactory bulb dataset       | Stereo-seq          | 1          | 10000           | 26145         | [link](https://doi.org/10.1016/j.cell.2022.04.003)                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| mouse brain cerebellum dataset     | Slide-seq           | 1          | 25551           | 20141         | [link](https://www.science.org/doi/10.1126/science.aaw1219)                                                                                                                                                                                                                                                                                                                                                                                                                                      |



- DISCLAIMER: mouse olfactory bulb dataset, mouse kidney coronal dataset, mouse brain sagittal dataset, mouse somatosensory cortex dataset, mouse olfactory bulb dataset and mouse brain cerebellum dataset are converted from [**this repository**](https://github.com/acheng416/Benchmark-CTCM-ST.git), which may have some mistakes. If you use those datasets mentioned above in your experiments, you should consider citing [**this paper**](https://academic.oup.com/bib/article/doi/10.1093/bib/bbac475/6835380).


# guide
Clone this repository into your project root first.
Then you can simply using all the dataset by the code below.

``` python
from st_datasets.dataset import get_data, dataset_you_need  # or you can use * to import all

adata, num = get_data(dataset_func=dataset_you_need, **dataset_specific_args)
```
`get_data()` will return a `sc.Anndata` type st data and the number of clusters in the dataset. 

# more datasets
squidpy provide some [annotated ST datasets](https://squidpy.readthedocs.io/en/stable/api.html#module-squidpy.datasets).  
You can use those datasets by the following code.

``` python
import squidpy as sq

adata = sq.dataset.dataset_you_need()
```
It may take some time to download the datasets when you first use them.

# some tools
- build graph: *(comming soon)*  
You can build adjacency matrix by using `utils.build_graph()`.  
We provide two methods to build the graph, KNN or build by radius.  
We support to return pyg `Data` format data with edge_index and x attribute or `np.array` format adjacency matrix.  

- show data distribution:  
You can draw a data distribution map by using `utils.show_distrib_map()`.  
It can draw a histgram or a violin plot using `seaborn`.  
You can customize your plot by adding seabon parameters.
