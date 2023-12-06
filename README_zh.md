# st_datasets 

[English](./README.md) | 中文简体

<details>
<summary>数据集详细信息</summary>
<table>
    <thead>
        <tr>
            <th>名称</th>
            <th>测序技术</th>
            <th>切片数量</th>
            <th>spots/细胞 个数</th>
            <th>基因个数</th>
            <th>来源</th>
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

\* 这部分数据集由[**this repository**](https://github.com/acheng416/Benchmark-CTCM-ST.git)数据转化打包为h5ad数据，如果您的实验使用到了上述数据，请考虑引用[**this paper**](https://academic.oup.com/bib/article/doi/10.1093/bib/bbac475/6835380).

\^ 由Stereo-seq测序获得的Mouse olfactory bulb dataset在新提交中被更新，*新版本数据**不提供**标签*，但原始数据依然保存在[**the huggingface repo**](https://huggingface.co/datasets/han-shu/st_datasets)中，如果您需要使用原数据，请使用先前提交版本(commit [**#c56d877**](https://github.com/hannshu/st_datasets/tree/c56d877a001071cb7b2cb4c222491ce20e026c22)
)的`st_datasets`!

</details>

## 使用说明
首先将本仓库克隆到您项目的根目录下，加下来使用一下的代码即可使用数据。

``` python
import st_datasets as stds

# 查找您需要使用的数据
stds.datasets()

# 读取并使用数据集
adata, num = stds.get_data(dataset_func=stds.dataset_you_need, **dataset_specific_args)
```
`get_data()` 将返回一个 `sc.Anndata` 类型的空间转录组数据以及一个该数据集的聚类类别个数。  
- `st_datasets`将会在第一次启动时生成一个代理配置文件。如果您不需要使用代理下载数据集，请再次
运行您的脚本，如果您需要配置代理，请依照指示设置代理。若您无法下载数据集，我们还将数据存放在[**Baidu Netdisk**](https://pan.baidu.com/s/1eMVnLnJvx17Q8NmGgikuZA?pwd=k3k5)上，请在下载完成后将`data`文件夹
放置在用户根目录下。

### 更多数据集
squidpy 也提供了一些数据集 [annotated ST datasets](https://squidpy.readthedocs.io/en/stable/api.html#module-squidpy.datasets).  
您可以通过以下代码获取该数据集中的数据。

``` python
import squidpy as sq

adata = sq.dataset.dataset_you_need()
```
在第一次使用时将自动下载数据集。

## 工具
### 预处理 
`st_datasets` 目前支持以下预处理方法:  
- `stds.pp.build_graph()` 生成图。   
- `stds.pp.convert_edge_to_adj()` and `stds.pp.convert_adj_to_edge()` 允许您在边集合和邻接矩阵之间进行转换。   
- `stds.pp.get_not_adjacency_pair()` 获取不邻接的点对。   
- `stds.pp.add_self_loop()` and `stds.pp.delete_self_loop()` 为图添加/删除自环的边。   
- `stds.pp.concat_adjacency_matrix()` 集合多张图。  
- `stds.pp.conv_to_one_hot()` 将标签转化为独热码。  
- `stds.pp.k_hop_adj()` 生成k跳邻接矩阵。 

### 绘图
`st_datasets` 目前支持以下绘图方法:  
- `stds.pl.show_distrib_map()` 绘制分布图。  
- `stds.pl.visualize_graph()` 绘制空间分布图。  

### 聚类
您可以使用 `stds.cl.evaluate_embedding()` 进行聚类  
`st_datasets` 目前支持以下聚类方法:  
- R package mclust-EEE algorithm  
- K-Means clustering method
- Birch clustering method  
- Gaussian Mixture clustering method  
- Agglomerative clustering method  


`st_datasets` 目前支持以下聚类评价指标:  
- ARI
- NMI
- AMI
- V measure score
- silhouette score (不需要标签)
