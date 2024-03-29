import os
import gc
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def start() -> None:
    from globals import H5AD_DIR, Clustresults
    
    d = 'Immune'
    # Load batch correct data 
    dest = f"{H5AD_DIR}/adata_final_{d}_cca_features.h5ad"
    if os.path.exists(dest):  ##checks if the batch-corrected data already exists
        print("Load batch corrected data...")
        print(dest)
        adata = sc.read_h5ad(dest) 
    else:
        return 
    
    resolutions = ['test_leiden_n15_r1.0', 'test_leiden_n15_r2.0']
    
    # UMAP 
    for r in resolutions:
        fig, ax = plt.subplots(figsize=(7, 7))
        adata.obsm["X_umap"] = adata.obsm["X_umap_neighbors_n15"]
        sc.pl.umap(adata,
                   use_raw=False,          # gene expression
                   color=r,
                   wspace=0.3,
                   #ncols=3,
                   neighbors_key='neighbors_15',
                   show=False,
                   ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title(f'Astrocytes Umap Resolution = {r[-3:0]} (Leiden)')
        fig.savefig(fname=f"{Clustresults}/{d}_umap_{r}.png",
                    dpi=300,
                    bbox_inches='tight')
        plt.close(fig)


start()     
