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
    from globals import H5AD_DIR, CLUSTER_FUSION_DIR, Clustresults, IMMUNE_MARKERS_CSV
    from globals import DATASETS, lineage_resolution_tests, ENV
    from helper_functions import random_colors
    
    for d in DATASETS:
        # The code iterates over datasets stored in the datasets_divided list
        # Load cluster fusion file. The first columns contain higher resolution culsters
        # The fusion will be done from bottom to top of the cluster tree.
        # Thus, if a cell is to belong to a cluster from a high res, it will not be part of a cluster form a lower res.
        
        # TODO: Adapt to import . tsv of known biomarkers
        dest = f"{ENV}/{d}_Markers.csv"
        if os.path.exists(dest):
            print("Load cluster fusion file...")
            print(dest)
            cluster_df = pd.read_csv(dest,
                                     sep='\t',
                                     header=0,
                                     dtype='category',
                                     index_col=None,
                                     na_filter=False)
            print(cluster_df)
        else:
            continue
        
        # STEPS TODO
        # 1. Load feature data (clusters). Only contains HVG genes (4000 genes).
        # 2. Load entire data (with all genes). File 'raw_norm_annot.h5ad'
        # 3. Integrate metada
        # 4. Do DGE with Wilcoxon 
        
        # 1
        # Load batch correct data
        dest = f"{H5AD_DIR}/adata_final_{d}_cca_features.h5ad"
        if os.path.exists(dest):    # checks if the batch-corrected data already exists
            print("Load batch corrected data...")
            print(dest)
            adata = sc.read_h5ad(dest)
            
            print(adata)
            
        else:
            continue
        
        # 2
        # TODO
        # Load 'raw_norm_annot.h5ad'
        
        for u in DATASETS:
            data = sc.read_h5ad(f"adata_final_{u}_raw_norm_annot.h5ad")
            print(data)
            

        # 3
        # Refator code
        # Transfer over the metadata
        data['raw_norm_annot'].obsm['X_umap'] = adata.obsm['X_umap'].copy()
        data['adata_raw_norm'].obsm['X_tsne'] = adata.obsm['X_tsne'].copy()
        
        resolution = lineage_resolution_tests[d]
        i = '15'
        for r in resolution:
            if f'X_umap_neighbors_n{i}' in data[d].obsm.keys():
                data['adata_raw_norm'].obsm[f'X_umap_neighbors_n{i}'] = data[d].obsm[f'X_umap_neighbors_n{i}'].copy()
            if f'test_leiden_n{i}_r{r}' in data[d].obs.keys():
                data['adata_raw_norm'].obs[f'test_leiden_n{i}_r{r}'] = data[d].obs[f'test_leiden_n{i}_r{r}'].copy()
            if f'dendrogram_leiden_n{i}_r{r}' in data[d].uns.keys():
                data['adata_raw_norm'].uns[f'dendrogram_leiden_n{i}_r{r}'] = data[d].uns[f'dendrogram_leiden_n{i}_r{r}'].copy()
        
        if 'leiden_fusion' in data[d].obs.keys():
            data['adata_raw_norm'].obs['leiden_fusion'] = data[d].obs['leiden_fusion'].copy()
        if 'dendrogram_leiden_fusion' in data[d].uns.keys():
            data['adata_raw_norm'].uns['dendrogram_leiden_fusion'] = data[d].uns['dendrogram_leiden_fusion'].copy()

        # Free some memory
        gc.collect()
        
        # 4
        # Find marker genes
        print("Find marker genes ...")
        # Marker genes for cluster fusions
        sc.tl.rank_genes_groups(adata=data['adata_raw_norm'],
                                n_genes=None,
                                groupby='leiden_fusion',
                                method='wilcoxon',
                                use_raw=False,
                                key_added='rank_genes_groups_leiden_fusion',
                                pts=True)
        
        # Free some memory
        gc.collect()
        
        # Save data
        print("Save gene rank data...")
        dest = f"{H5AD_DIR}/adata_final_{d}_raw_norm_ranked.h5ad"
        print(dest)
        print(data['adata_raw_norm'])
        data['adata_raw_norm'].write_h5ad(dest, compression='gzip')

       


start()
