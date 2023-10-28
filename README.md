# st_datasets 


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
            <th>mouse olfactory bulb dataset*</th>
            <td>Stereo-seq</td>
            <td>1</td>
            <td>10000</td>
            <td>26145</td>
            <td>converted</td>
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

### more datasets
squidpy provide some [annotated ST datasets](https://squidpy.readthedocs.io/en/stable/api.html#module-squidpy.datasets).  
You can use those datasets by the following code.

``` python
import squidpy as sq

adata = sq.dataset.dataset_you_need()
```
It may take some time to download the datasets when you first use them.

## tools
- build graph: *(comming soon)*  
You can build adjacency matrix by using `utils.build_graph()`.  
We provide two methods to build the graph, KNN or build by radius.  
We support to return pyg `Data` format data with edge_index and x attribute or `np.array` format adjacency matrix.  

- show data distribution:  
You can draw a data distribution map by using `utils.show_distrib_map()`.  
It can draw a histgram or a violin plot using `seaborn`.  
You can customize your plot by adding seabon parameters.
