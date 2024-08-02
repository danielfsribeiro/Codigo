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
        # The code iterates over datasets stored in the datasets_divided list
        # Load cluster fusion file. The first columns contain higher resolution culsters
        # The fusion will be done from bottom to top of the cluster tree.
        # Thus, if a cell is to belong to a cluster from a high res, it will not be part of a cluster form a lower res.
        
        
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
        adata_raw = sc.read_h5ad(f"{H5AD_DIR}/adata_final_{d}_raw_norm_annot.h5ad")
        print(adata_raw)
        

        # 3
        # Refator code
        # Transfer over the metadata
        adata_raw.obsm['X_umap'] = adata.obsm['X_umap'].copy()
        adata_raw.obsm['X_tsne'] = adata.obsm['X_tsne'].copy()
        
        resolution = lineage_resolution_tests[d]
        i = '15'
        if f'X_umap_neighbors_n{i}' in adata.obsm.keys():
            adata_raw.obsm[f'X_umap_neighbors_n{i}'] = adata.obsm[f'X_umap_neighbors_n{i}'].copy()

        if f"neighbors_{i}" in adata.obs.keys():
            adata_raw.obs[f'neighbors_{i}'] = adata.obs[f'X_umap_neighbors_n{i}'].copy()


        
        for r in resolution:
            
            #for item in data[d].obs.keys():
            #    'leiden' in item:
            #        # DO copy

            if f'test_leiden_n{i}_r{r}' in adata.obs.keys():
                adata_raw.obs[f'test_leiden_n{i}_r{r}'] = adata.obs[f'test_leiden_n{i}_r{r}'].copy()
            
            if f'leiden_n{i}_r{r}' in adata.obs.keys():
                adata_raw.obs[f'leiden_n{i}_r{r}'] = adata.obs[f'leiden_n{i}_r{r}'].copy()

            if f'dendrogram_leiden_n{i}_r{r}' in adata.uns.keys():
                adata_raw.uns[f'dendrogram_leiden_n{i}_r{r}'] = adata.uns[f'dendrogram_leiden_n{i}_r{r}'].copy()
        
        # TODO
        ### Tens de remover o comentátio destas linhas seguintes, as linhas 91-95.
        ### No ficheiro anterior "Cluster_fusion.py" adicionas .obs['leiden_fusion'] ao ficheiro 'features' nas linhas 83-84.
        ### Depois tens de copiar toda essa info para o ficherio 'raw', que é feito nessas linhas seguintes.
        ### Penso que depois já tudo funcione.
        
        if 'leiden_fusion' in adata_raw.obs.keys():
            adata_raw.obs['leiden_fusion'] = adata.obs['leiden_fusion'].copy()

        if 'dendrogram_leiden_fusion' in adata_raw.uns.keys():
            adata_raw['adata_raw_norm'].uns['dendrogram_leiden_fusion'] = adata.uns['dendrogram_leiden_fusion'].copy()

        # Free some memory
        gc.collect()
        
        # 4
        #log transformation of the data
        sc.pp.log1p(adata_raw)

        # Find marker genes
        print("Find marker genes ...")
        # Marker genes for cluster fusions
        sc.tl.rank_genes_groups(adata=adata_raw,
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
        adata_raw.write_h5ad(dest, compression='gzip')

       


start()
