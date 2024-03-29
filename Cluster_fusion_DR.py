import os
import gc
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

### Defines
#Good_clusters_r2 = ['1', '10', '15', '24', '30']      # .obs['test_leiden_n15_r2']
#Good_clusters_r1 = ["6", "7", "9"]
#Bad_clusters_r1 = ["0", "1", "2", "3", "4,", "5", "8"]


# Convert old cluter nr into new cluster nr
dict_good_cells = {
    '1': '2.1',
    '10': '2.10',
    "15": "2.15",
    "24": "2.24",
    "30": "2.30"
}
dict_bad_cells = {
    '0': '1.0',
    '1': '1.1',
    '2': '1.2',
    '3': '1.3',
    '4': '1.4',
    '5': '1.5',
    '8': '1.8',
}

### Load data
# Use adata_cca_features.h5ad --> must contain info about r1 and r2

# New code to refactor
small_names = {
    'Astrocyte': 'Ast',
    'Immune': 'Imm',
}


def start() -> None:
    # Import working directories
    from globals import H5AD_DIR, CLUSTER_FUSION_DIR
    from globals import DATASETS
    
    for d in DATASETS:   
        #The code iterates over datasets stored in the datasets_divided list
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


        #TODO
        
        clustered_cells = []                 # Save cells that were already annotated in a higher resolution
        print(adata)
        for res in cluster_df.columns:
            key = f'test_leiden_n15_{res}.00'      # recupera a resolução que precisamos olhar
            categories = [str(s) for s in cluster_df.loc[:, res].cat.categories.to_list()]
            print(f"Clusters to include res={key}: {' '.join(categories)}")     #Debug
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
            
        
        #TODO
        # Cells not clustered - label them as cluster 'NA'
        unclustered_cells = adata[~adata.obs_names.isin(clustered_cells)].obs_names.to_list()
        print(f"Number of cells not selected: {len(unclustered_cells)}")
        # Update AnnData
        adata.obs.loc[unclustered_cells, 'leiden_fusion'] = f"{small_names[d]}.NA"
        adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype(dtype='category')
        # Check if any new cluster has less than 10 cells. If so, label them as 'NA'
        clusters = adata.obs['leiden_fusion'].cat.categories.to_list()
        clusters.sort()
        for c in clusters:
            mask = adata.obs['leiden_fusion'] == c
            if mask.sum() < 10:
                print(f"Cluster {c} has less than 10 cells! Labeling cells as 'NA'")
                adata.obs.loc[mask, 'leiden_fusion'] = f"{small_names[d]}.NA"
        # Reset clusters
        adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype(dtype='category')
        adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].cat.remove_unused_categories()

        
        ## Convert cluster names

        #Células pertencentes aos clusters com resolução final aceitável de 2.00 (Goodcells)
        #for key, value in dict_good_cells.items():
        #    mask_3 = adata.obs['leiden_fusion'].isin([key])
        #    adata.obs.loc[mask_3, 'leiden_fusion'] = value

        #Células pertencentes aos clusters com resolução final aceitável de 1.00 (Badcells)
        #for key, value in dict_bad_cells.items():
        #    mask_4 = adata.obs['leiden_fusion'].isin([key])
        #    adata.obs.loc[mask_4, 'leiden_fusion'] = value
        
        #columns_to_keep = ["Fusion_cluster"] #apagar obs: "leiden_n15", "test_leiden_n15"
        #adata.obs = adata.obs[columns_to_keep]
        #adata.uns.clear()
        
        #print("Goodcells:")
        #print(Goodcells.obs["test_leiden_n15_r2.00"].to_list())
        #print("Badcells:")
        #print(Goodcells.obs["test_leiden_n15_r1.00"].to_list())
        #print("adata")
        #print(adata.obs["Fusion_cluster"].to_list())
        #print(adata.obs)
        #print(adata.uns)
        
        # UMAP 
        #fig = plt.figure(figsize=(7, 7), dpi=300)
        fig, ax = plt.subplots(figsize=(7, 7))
        
        adata.obsm["X_umap"] = adata.obsm["X_umap_neighbors_n15"]
        sc.pl.umap(adata,
                   use_raw=False,          # gene expression
                   color='test_leiden_n15_r1.00',
                   wspace=0.3,
                   #ncols=3,
                   neighbors_key='neighbors_15',
                   show=False,
                   ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title('Immune Umap Resolution = 1.00 (Leiden)')
        fig.savefig(fname="/home/nolanreardson/Bioinformatics_Env_Windows/Astrocytes_imune_Cells/Clustresults/Immune_umap_test_leiden_n15_r1.00.png",
                    dpi=300,
                    bbox_inches='tight')
        plt.close(fig)

        #fig = plt.figure(figsize=(7, 7), dpi=300)
        fig, ax = plt.subplots(figsize=(7, 7))

        sc.pl.umap(adata,
                   use_raw=False,          # gene expression
                   color='test_leiden_n15_r2.00',
                   wspace=0.3,
                   #ncols=3,
                   neighbors_key='neighbors_15',
                   show=False,
                   ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title('Immune Umap Resolution = 2.00 (Leiden)')
        fig.savefig(fname="/home/nolanreardson/Bioinformatics_Env_Windows/Astrocytes_imune_Cells/Clustresults/Immune_umap_test_leiden_n15_r2.00.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

        #fig = plt.figure(figsize=(7, 7), dpi=300)
        fig, ax = plt.subplots(figsize=(7, 7))
        
        sc.pl.umap(adata,
                   use_raw=False,          # gene expression
                   color='leiden_fusion',
                   wspace=0.3,
                   #ncols=3,
                   neighbors_key='neighbors_15',
                   show=False,
                   ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title('Immune Umap Fusion Cluster (Leiden)')
        fig.savefig(fname="/home/nolanreardson/Bioinformatics_Env_Windows/Astrocytes_imune_Cells/Clustresults/Immune_umap_leiden_fusion.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

        
        #group_1 = adata.obs["Fusion_cluster"].str.startswith("1.")
        #group_2 = adata.obs["Fusion_cluster"].str.startswith("2.")
        
        # Perform UMAP for each group separately
        #sc.pp.neighbors(group_1)
        #sc.tl.umap(group_1)

        #sc.pp.neighbors(group_2)
        #sc.tl.umap(group_2)
        
        # Plot UMAPs for both groups
        #sc.pl.umap(group_1)
        #sc.pl.umap(group_2)
 
 ## Convert fusion cluster to categories: 
 
 ### How to convert to categories:
# data[d].obs[f'{es}_rank'] = None
# code... code ...
#   data[d].obs.loc[top_cells, f'{es}_rank'] = 'top'
#   data[d].obs.loc[low_cells, f'{es}_rank'] = 'low'
# 
# When column is finished: 
# data[d].obs[f'{es}_rank'] = data[d].obs[f'{es}_rank'].astype(dtype='category')
#adata.obs['fusion_cluster'] = adata.obs['fusion_cluster'].astype("category")
      
start()     
        
        # ....
        
        # Delete other clusters (manter apenas "Fusion cluster", eliminar "...leiden_n15..." das secções .obs e .uns )
        # Salvar como ... outro nome --> sc.write_h5ad(adata, output_filename) print(f"Altered data saved to {output_filename}")
        
        

# ...
# TODO 3: print of categories names
#print(small_adata.obs['fusion_cluster'].cat.categories)

### Output
# Fazer um output de umaps da resolução test_leiden_r1 , r2, r1_2 (fusion of r1 and r2) 
#print(small_adata.obs['test_leiden_n15_r1'].cat.categories) 


          

