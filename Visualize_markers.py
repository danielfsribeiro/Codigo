import os
import gc
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# New code to refactor
small_names = {
    'Astrocyte': 'Ast',
    'Immune': 'Imm',
}


def start() -> None:
    # Import working directories
    from globals import H5AD_DIR, CLUSTER_FUSION_DIR, Clustresults
    from globals import DATASETS, lineage_resolution_tests, ENV
    from helper_functions import random_colors
    
    for d in DATASETS:
        
        
        # STEPS TODO
        # 1. Load literature marker files.
        dest = f"{ENV}/{d}_Markers.csv"
        if os.path.exists(dest):
            print("Load cluster fusion file...")
            print(dest)
            cluster_df = pd.read_csv(dest,
                                     sep=';',
                                     header=0,
                                     dtype='category',
                                     index_col=None,
                                     na_filter=False)
            print(cluster_df)
        else:
            continue

        # 2. Load entire data (with all genes). File 'raw_norm_annot.h5ad'
        dest = f"{H5AD_DIR}/adata_final_{d}_raw_norm_annot.h5ad"
        if os.path.exists(dest):    # checks if the batch-corrected data already exists
            print("Load batch corrected data...")
            print(dest)
            adata_raw = sc.read_h5ad(dest)
            
            print(adata_raw)
                        
        else:
            continue

        # 3. Visualize DGE markers (sc.pl.(...dotplot rank...)) 
        inj = adata_raw.obs["injury"]
        

        
        sc.pl.dotplot(adata = adata_raw, var_names= cluster_df, groupby = inj , use_raw=False, log=False, dendrogram=True)

        # 4. Vizualize literature markers (sc.pl.dotplot) -> dar a lista de genes de marcadores e a lista de tipos de celulas

        cell_types = {}
        cell_types['M1 macrophage'] = ['CD80', 'CD86']
        cell_types['M2 macrophage'] = ['MRC1', 'IL-10']


        # 5. Vizualize literature markers (pl.umap), do a umap for each gene
        

        

       
       

      
