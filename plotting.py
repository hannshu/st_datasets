import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def show_distrib_map(adata, title='data_distribution_map', format='hist', **args):
    assert format in ['hist', 'violin'], \
        ">>> ERROR: This return data type is not provided now."

    if ('highly_variable' in adata.var):
        adata = adata[:, adata.var['highly_variable']]

    val = adata.X.todense().reshape(adata.X.shape[0]*adata.X.shape[1]).tolist()[0]
    if ('violin' == format):
        df = pd.DataFrame()
        df['gen_exp'] = val
        df['name'] = [title]*len(val)
        sns.violinplot(data=df, **args) 
        plt.title(title)
    else:
        sns.histplot(val, **args)
        plt.title(title)
    plt.savefig(title+'.png', dpi=600)