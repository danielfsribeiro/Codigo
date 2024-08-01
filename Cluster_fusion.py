import os
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
    from globals import DATASETS
    from helper_functions import random_colors
    
    for d in DATASETS:
        # The code iterates over datasets stored in the datasets_divided list
        # Load cluster fusion file. The first columns contain higher resolution culsters
        # The fusion will be done from bottom to top of the cluster tree.
        # Thus, if a cell is to belong to a cluster from a high res, it will not be part of a cluster form a lower res.
        dest = f"{CLUSTER_FUSION_DIR}/{d}_fusion.tsv"
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
        
        # Load batch correct data
        dest = f"{H5AD_DIR}/adata_final_{d}_cca_features.h5ad"
        if os.path.exists(dest):    # checks if the batch-corrected data already exists
            print("Load batch corrected data...")
            print(dest)
            adata = sc.read_h5ad(dest)
        else:
            continue

        clustered_cells = []                 # Save cells that were already annotated in a higher resolution
        print(adata)
        for res in cluster_df.columns:
            key = f'test_leiden_n15_{res}'      # recupera a resolução que precisamos olhar
            categories = [str(s) for s in cluster_df.loc[:, res].cat.categories.to_list()]
            print(f"Clusters to include res={key}: {' '.join(categories)}")     # Debug
            # Which cells belonging to selected clusters
            temp = adata.obs.loc[:, key].copy()       # Para cada celula a qual cluster pertence na resolução leidein_rxx
            mask = temp.isin(cluster_df.loc[:, res])    # Dentro da resolução, quais pertencem aos clusters da coluna em analise em cluster_df
            temp = adata[mask].copy()                 # Just to be safe
            # Remove cells from the already clustered list
            mask2 = temp.obs_names.isin(clustered_cells)    # Ver se as células já foram clustered before
            temp = temp[~mask2].copy()                      # Se sim, remove-las da anotação de 'leidein_fusion'
            print(f"Number of cells selected: {temp.shape[0]}")
            clustered_cells.extend(temp.obs_names)  # Add new cells to list of already clustered cells
            # Update AnnData
            adata.obs.loc[temp.obs_names, 'leiden_fusion'] = [f"{small_names[d]}.{str(res[1:])}.{c}" for c in temp.obs.loc[:, key].to_list()]
        
        
        # Cells not clustered - label them as cluster 'NA'
        unclustered_cells = adata[~adata.obs_names.isin(clustered_cells)].obs_names.to_list()
        print(f"Number of cells not selected: {len(unclustered_cells)}")
        print(adata[~adata.obs_names.isin(clustered_cells)].obs.loc[:, ['test_leiden_n15_r2.0', 'test_leiden_n15_r1.0']])
        # Update AnnData
        adata.obs.loc[unclustered_cells, 'leiden_fusion'] = f"{small_names[d]}.NA"
        adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype(dtype='category')
        # Check if any new cluster has less than 10 cells. If so, label them as 'NA'
        clusters = adata.obs['leiden_fusion'].cat.categories.to_list()
        clusters.sort()
        print(clusters)
        for c in clusters:
            mask = adata.obs['leiden_fusion'] == c
            if mask.sum() < 10:
                print(f"Cluster {c} has less than 10 cells! Labeling cells as 'NA'")
                adata.obs.loc[mask, 'leiden_fusion'] = f"{small_names[d]}.NA"
        # Reset clusters
        adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype(dtype='category')
        adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].cat.remove_unused_categories()
        
        # UMAP
        fig, ax = plt.subplots(figsize=(7, 7))
        
        colors = random_colors(len(adata.obs['leiden_fusion'].cat.categories.to_list()))
        sc.pl.umap(adata,
                   use_raw=False,          # gene expression
                   color='leiden_fusion',
                   wspace=0.3,
                   #ncols=3,
                   palette=colors,
                   neighbors_key='neighbors_15',
                   show=False,
                   ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title(f'{d} Umap Fusion Cluster (Leiden)')
        fig.savefig(fname=f"{Clustresults}/{d}_umap_leiden_fusion.png", dpi=300, bbox_inches='tight')
        plt.close(fig)


    # Save cca_features.h5ad
    # Save data - TODO
    dest = f"{H5AD_DIR}/adata_final_{d}_cca_features.h5ad"
    print(dest)
    adata.write_h5ad(dest, compression='gzip')


start()
