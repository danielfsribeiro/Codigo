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
        
        # TODO: Adapt to import . tsv of known biomarkers
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
        # 3. Visualize DGE markers (sc.pl.(...dotplot rank...)) 
        # 4. Vizualize literature markers (sc.pl.dotplot) -> dar a lista de genes de marcadores e a lista de tipos de celulas
        # 5. Vizualize literature markers (pl.umap), do a umap for each gene
        

        cell_types = {}
        cell_types['M1 macrophage'] = ['CD80', 'CD86']
        cell_types['M2 macrophage'] = ['MRC1', 'IL-10']

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
        adata_raw = sc.read_h5ad(f"{H5AD_DIR}/adata_final_{d}_raw_norm_annot.h5ad")
        print(adata_raw)
        

        # 3
        # Refator code
        # Transfer over the metadata
        adata_raw.obsm['X_umap'] = adata.obsm['X_umap'].copy()
        adata_raw.obsm['X_tsne'] = adata.obsm['X_tsne'].copy()
        
        resolution = lineage_resolution_tests[d]
        i = '15'
        if f'X_umap_neighbors_n{i}' in data[d].obsm.keys():
             data['adata_raw_norm'].obsm[f'X_umap_neighbors_n{i}'] = data[d].obsm[f'X_umap_neighbors_n{i}'].copy()
    
        for r in resolution:
            
            #for item in data[d].obs.keys():
            #    'leiden' in item:
            #        # DO copy

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
