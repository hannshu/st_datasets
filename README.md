# st_datasets 

English | [中文简体](./README_zh.md)

<details>
<summary>Dataset details</summary>
<table>
    <thead>
        <tr>
            <th>dataset</th>
            <th>technology</th>
            <th>slices</th>
            <th>spots/cells</th>
            <th>genes</th>
            <th>source</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <th>DLPFC dataset</th>
            <td>10x Genomics Visium</td>
            <td>12</td>
            <td>3460 - 4789</td>
            <td>33538</td>
            <td>10x visium database</td>
        </tr>
        <tr>
            <th>human breast cancer dataset</th>
            <td>10x Genomics Visium</td>
            <td>1</td>
            <td>3798</td>
            <td>36601</td>
            <td>10x visium database</td>
        </tr>
        <tr>
            <th>mouse olfactory bulb dataset</th>
            <td>ST</td>
            <td>12</td>
            <td>231 - 282</td>
            <td>15284 - 16675</td>
            <td>stomicsDB</td>
        </tr>
        <tr>
            <th>mouse kidney coronal dataset*</th>
            <td>10x Genomics Visium</td>
            <td>1</td>
            <td>1438</td>
            <td>31053</td>
            <td>converted</td>
        </tr>
        <tr>
            <th>mouse brain sagittal dataset</th>
            <td>10x Genomics Visium</td>
            <td>4</td>
            <td>2696 - 3353</td>
            <td>31053</td>
            <td>stomicsDB</td>
        </tr>
        <tr>
            <th>mouse somatosensory cortex dataset*</th>
            <td>osmFISH</td>
            <td>1</td>
            <td>5328</td>
            <td>33</td>
            <td>converted</td>
        </tr>
        <tr>
            <th>mouse olfactory bulb dataset^</th>
            <td>Stereo-seq</td>
            <td>1</td>
            <td>19109</td>
            <td>27106</td>
            <td>https://doi.org/10.1016/j.cell.2022.04.003</td>
        </tr>
        <tr>
            <th>mouse brain cerebellum dataset*</th>
            <td>Slide-seq</td>
            <td>1</td>
            <td>25551</td>
            <td>20141</td>
            <td>converted</td>
        </tr>
        <tr>
            <th>Mouse Organogenesis Spatiotemporal Transcriptomic Atlas (E9.5)</th>
            <td>Stereo-seq</td>
            <td>5</td>
            <td>4356 - 5931</td>
            <td>24238</td>
            <td>stomicsDB</td>
        </tr>
        <tr>
            <th>Zebrafish Embryogenesis Spatiotemporal Transcriptomic Atlas</th>
            <td>Stereo-seq</td>
            <td>1 (with 6 sections)</td>
            <td>13166</td>
            <td>26628</td>
            <td>stomicsDB</td>
        </tr>
    </tbody>
</table>

\* Dataset are converted from [**this repository**](https://github.com/acheng416/Benchmark-CTCM-ST.git). If you use those datasets mentioned above in your experiments, you should consider citing [**this paper**](https://academic.oup.com/bib/article/doi/10.1093/bib/bbac475/6835380).

\^ Mouse olfactory bulb dataset obtained by Stereo-seq has been updated in the new commit, *the new version data **do not** have annotation*, but the original data is not removed from [**the huggingface repo**](https://huggingface.co/datasets/han-shu/st_datasets), if you need the old version data, please use the former commit (commit [**#c56d877**](https://github.com/hannshu/st_datasets/tree/c56d877a001071cb7b2cb4c222491ce20e026c22)
) of `st_datasets`!

</details>

## guide
Clone this repository into your project root first.
Then you can simply using all the dataset by the code below.

``` python
import st_datasets as stds

# find datasets you need
stds.datasets()

# load and use dataset
adata, num = stds.get_data(dataset_func=stds.dataset_you_need, **dataset_specific_args)
```
`get_data()` will return a `sc.Anndata` type st data and the number of clusters in the dataset.  
- Please note that `st_datasets` may generate a proxy file the first time it is set up, which
may help people who cannot directly download datasets. 
If you do not need a proxy, just rerun your script again, and the error will disappear.
For those who cannot access the remote dataset and cannot set a proxy, we provide [**Baidu Netdisk**](https://pan.baidu.com/s/1eMVnLnJvx17Q8NmGgikuZA?pwd=k3k5) 
to download the dataset. Please place the `data` folder in your user root directory.

- We encourage you to find the data you need from [**stomicsdb**](https://db.cngb.org/stomics/) or other st database. You can easily download data
from an outside database by `stds.get_data(stds.get_outside_data, url=data_url)`. `st_datasets` will download your data to `~/data/user_download/`.
You can access your data again by providing the data_url.

### more datasets
squidpy provide some [annotated ST datasets](https://squidpy.readthedocs.io/en/stable/api.html#module-squidpy.datasets).  
You can use those datasets by the following code.

``` python
import squidpy as sq

adata = sq.dataset.dataset_you_need()
```
It may take some time to download the datasets when you first use them.

## tools
### Preprocessing 
You can use `st_datasets` to do preprocessing such as:  
- `stds.pp.build_graph()` build spatial graph.  
- `stds.pp.convert_edge_to_adj()` and `stds.pp.convert_adj_to_edge()` can convert edge list to adjacency matrix, verse vise.  
- `stds.pp.get_not_adjacency_pair()` get non adjacency spot pair (not linked edges).  
- `stds.pp.add_self_loop()` and `stds.pp.delete_self_loop()` add or delete self loop of the graph.  
- `stds.pp.concat_adjacency_matrix()` concate multiple adjacency matrix.  
- `stds.pp.conv_to_one_hot()` convert category label to one hot.  
- `stds.pp.k_hop_adj()` generate k hop adjacency matrix. 

### Plotting
You can use `st_datasets` to do plotting such as:  
- `stds.pl.show_distrib_map()` draw distribution image.  
- `stds.pl.visualize_graph()` draw spatial graph.  

### Clustering
You can use `st_datasets` to do clustering by `stds.cl.evaluate_embedding()`.  
We currently support clustering methods such as:  
- R package mclust-EEE algorithm  
- K-Means clustering method
- Birch clustering method  
- Gaussian Mixture clustering method  
- Agglomerative clustering method  


We also support evaluation methods such as:  
- ARI
- NMI
- AMI
- V measure score
- silhouette score (unsupervised)
