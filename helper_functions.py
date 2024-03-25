# # Helper functions

from typing import Union
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.sparse import issparse


def plot_umap(data: sc.AnnData,
              dest: str,
              color: Union[str, list[str]],
              name_suffix: Union[str, None] = None,
              neighbors_key: Union[str, None] = None,
              multipanel: bool = False,
              multipanel_main_title: Union[str, None] = None,
              palette=None,
              **kwargs) -> None:
    ax = sc.pl.umap(data,
                    use_raw=False,
                    color=color,
                    palette=palette,
                    wspace=0.3,
                    ncols=3,
                    neighbors_key=neighbors_key,
                    show=False,
                    **kwargs)
    if multipanel:
        plt.suptitle(f"UMAP - {multipanel_main_title} - {name_suffix}")
        plt.savefig(dest + f"/UMAP - {multipanel_main_title} - {name_suffix}.pdf", bbox_inches="tight", dpi=600)
    else:
        # Dynamic detemrine legend columns
        numitems = len(list(ax.get_legend_handles_labels()[0]))
        nrows = 10
        ncols = int(np.ceil(numitems / float(nrows)))
        if ncols == 0:
            ncols = 1
        ax.legend(loc='center left', ncol=ncols, bbox_to_anchor=(1.05, 0.5), frameon=False)
        if 'projection' in kwargs.keys():
            if kwargs['projection'] == '3d':
                print('3d lines')
                ax.grid(True)
        plt.title(f"UMAP - {color} - {name_suffix}")
        plt.savefig(dest + f"/UMAP - {color} - {name_suffix}.pdf", bbox_inches="tight", dpi=600)
    plt.close()


def plot_tsne(data: sc.AnnData,
              dest: str,
              color: Union[str, list[str]],
              name_suffix: Union[str, None] = None,
              neighbors_key: Union[str, None] = None,
              multipanel: bool = False,
              multipanel_main_title: Union[str, None] = None,
              palette=None,
              **kwargs) -> None:
    ax = sc.pl.tsne(data,
                    use_raw=False,
                    color=color,
                    palette=palette,
                    alpha=0.5,
                    wspace=0.3,
                    ncols=3,
                    neighbors_key=neighbors_key,
                    show=False,
                    **kwargs)
    if multipanel:
        if name_suffix is None:
            plt.suptitle(f"tSNE - {multipanel_main_title}")
            plt.savefig(dest + f"/tSNE - {multipanel_main_title}.pdf", bbox_inches="tight")
        else:
            plt.suptitle(f"tSNE - {multipanel_main_title} - {name_suffix}")
            plt.savefig(dest + f"/tSNE - {multipanel_main_title} - {name_suffix}.pdf", bbox_inches="tight")
    else:
        # Dynamic detemrine legend columns
        numitems = len(list(ax.get_legend_handles_labels()[0]))
        nrows = 10
        ncols = int(np.ceil(numitems / float(nrows)))
        if ncols == 0:
            ncols = 1
        ax.legend(loc='center left', ncol=ncols, bbox_to_anchor=(1.05, 0.5), frameon=False)
        if name_suffix is None:
            plt.title(f"tSNE - {color}")
            plt.savefig(dest + f"tSNE - {color}.pdf", bbox_inches="tight")
        else:
            plt.title(f"tSNE - {color} - {name_suffix}")
            plt.savefig(dest + f"tSNE - {color} - {name_suffix}.pdf", bbox_inches="tight")
    plt.close()


def sparsity(x) -> float:
    if issparse(x):
        nz = x.count_nonzero()
        return 1.0 - float(nz) / float(np.prod(x.shape))
    else:
        nz = np.count_nonzero(x)
        return 1.0 - float(nz) / float(x.size)


def get_significant_genes(adata: sc.AnnData,                 # AnnData
                          key: str,                          # Key used to store the ranking results in adata.uns
                          n_genes: int = 5,                  # Number of genes to show
                          pval_cutoff: float = 0.05,         # pvals_adj cutoff to consider
                          expression_cutoff: float = 0.4,    # Expression cutoff to consider
                          #fc_cutoff: float = 1.0             # LogFC cutoff  (TOO RESTRICTIVE)
                          ) -> list:
    '''\
        Collect n_genes that have highest logFC and pvals_adj < 'pval_cutoff'
        and expression cutoff > 'expression_cutoff'
        key: `str`
            Key used to store the ranking results in adata.uns
        n_genes: `int` = 5
            Number of genes to show
        pval_cutoff: `float` = 0.05
            pvals_adj cutoff to consider
        expression_cutoff: `float` = 0.4
            Expression cutoff to consider (fraction of cell sexpressing the gene)
        fc_cutoff: `float` = 1.0
            Log FC cutoff to consider

        Returns
        -------
        list((genes, n genes added, bool for pval)). Tuple containing list of genes, the size and if they were selected from significant pval.
        Values will be sorted by highest logFC, then expression
    '''
    import pandas as pd
    # Get size col numbers
    n_col = len(pd.DataFrame(adata.uns[key]['names']).columns)
    gene_list = []
    for col in range(n_col):
        df = pd.DataFrame()
        df['names'] = pd.DataFrame(adata.uns[key]['names']).iloc[:, col]
        df['scores'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, col]
        df['logfoldchanges'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, col]
        df['pvals_adj'] = pd.DataFrame(adata.uns[key]['pvals_adj']).iloc[:, col]
        df['pts'] = pd.DataFrame(adata.uns[key]['pts']).loc[df['names'], str(col)].values
        df = df[df['pvals_adj'] < pval_cutoff].copy()
        df = df[df['pts'] >= expression_cutoff].copy()
        df.sort_values(by=['logfoldchanges', 'pts'], ascending=False, ignore_index=True, inplace=True)
        #df = df[df['logfoldchanges'] >= fc_cutoff].copy()
        values = df['names'].iloc[:n_genes].values
        gene_list.append([values.tolist(), len(values), True])  # Add True for significant pval
        # When no genes are found use highest score from within pts cutoff
        if gene_list[-1][1] == 0:
            df = pd.DataFrame()
            df['names'] = pd.DataFrame(adata.uns[key]['names']).iloc[:, col]
            df['scores'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, col]
            df['logfoldchanges'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, col]
            df['pvals_adj'] = pd.DataFrame(adata.uns[key]['pvals_adj']).iloc[:, col]
            df['pts'] = pd.DataFrame(adata.uns[key]['pts']).loc[df['names'], str(col)].values
            df = df[df['pts'] >= expression_cutoff].copy()
            df.sort_values(by=['logfoldchanges', 'pts'], ascending=False, ignore_index=True, inplace=True)
            values = df['names'].iloc[:n_genes].values
            gene_list[-1] = [values.tolist(), len(values), False]   # Add False for significant pval

    return gene_list
    
    
def random_colors(n: int, seed: Union[int, None] = 0) -> list:
    import numpy as np
    from scipy.interpolate import interp1d
    from matplotlib.colors import is_color_like
    
    # n must be smaller than 256 * 256 * 256
    if not n < 256 * 256 * 256:
        raise ValueError("n must be smaller than 256 * 256 * 256")
        
    colors = []
    if seed is None:
        randomizer = np.random.default_rng()
    else:
        randomizer = np.random.default_rng(seed=seed)
    scaler = interp1d([0, 255], [0.0, 1.0])
    i = 0
    while i < n:
        rgb = randomizer.choice(a=256, size=3)
        rgb = tuple(scaler(rgb).tolist())
        while (rgb in colors) or (not is_color_like(rgb)):
            rgb = randomizer.choice(a=256, size=3)
            rgb = tuple(scaler(rgb).tolist())
        colors.append(rgb)
        i += 1
    
    return colors
