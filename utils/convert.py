import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix


def bin2adata(file_path: str, distance: float = 50., loc_style: str = 'visium'):
    
    df = pd.read_csv(file_path, sep='\t')
    coors = df[['x', 'y']].to_numpy().astype('float32')

    # build coordination
    x_scale = [coors[:, 0].min(), coors[:, 0].max()]
    y_scale = [coors[:, 1].min(), coors[:, 1].max()]
    xs, ys = np.meshgrid(
        np.linspace(x_scale[0], x_scale[1], int((x_scale[1] - x_scale[0]) // distance)), 
        np.linspace(y_scale[0], y_scale[1], int((y_scale[1] - y_scale[0]) // distance))
    )
    if ('visium' == loc_style):
        xs[list(range(1, xs.shape[0], 2)), :] += distance / 2
    grid = np.array([xs.reshape(-1), ys.reshape(-1)]).T.astype('float32')
    
    print('Successfully loaded the bin file. Beginning to build the adata object; this process may take a considerable amount of time.')

    try:

        # use faiss to accelerate 
        import faiss

        index = faiss.IndexFlatL2(grid.shape[1])
        index.add(grid)

        _, df['nearest_idx'] = index.search(coors, 1) 

    except ImportError:

        # use sklearn as default
        from sklearn.neighbors import NearestNeighbors

        nbrs = NearestNeighbors(n_neighbors=1).fit(grid)
        _, df['nearest_idx'] = nbrs.kneighbors(coors)

    gene_df = pd.pivot_table(
        df, 
        values='MIDCounts', index='nearest_idx', columns='geneID', 
        aggfunc='sum', fill_value=0
    )

    # build adata
    adata = sc.AnnData(csr_matrix(np.float32(gene_df.to_numpy())))
    adata.var.index = gene_df.columns.to_list()
    adata.obs['total_count'] = gene_df.sum(axis=1)
    adata.obsm['spatial'] = grid[gene_df.index]

    print(f'>>> INFO: Successfully generate adata {adata.shape}')

    return adata