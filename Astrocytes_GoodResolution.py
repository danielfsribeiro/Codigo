import os
import gc
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

### Defines
 


# Convert old cluter nr into new cluster nr



### Load data
# Use adata_cca_features.h5ad --> must contain info about r1 and r2

Astrocytes_imune_Cells = "/mnt/c/Users/vasco/Env/Astrocytes_imune_Cells"

def start() -> None:
    import src.globals
    from src.globals import data, data2
  
    for d in ['Astrocyte']:    ##The code iterates over datasets stored in the datasets_divided list
        # Load batch correct data 
        dest = f"{Astrocytes_imune_Cells}/adata_final_{d}_cca_features.h5ad"
        if os.path.exists(dest):  ##checks if the batch-corrected data already exists
            print("Load batch corrected data...")
            print(dest)
            adata = sc.read_h5ad(dest)
                       
        else:
            continue  
        #test_leiden_n15_r0.60
        ## Convert cluster names
        
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
                   color='test_leiden_n15_r0.60',
                   wspace=0.3,
                   #ncols=3,
                   neighbors_key='neighbors_15',
                   show=False,
                   ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title('Astrocytes Umap Resolution = 0.60 (Leiden)')
        fig.savefig(fname="/home/nolanreardson/Bioinformatics_Env_Windows/Astrocytes_imune_Cells/Clustresults/Astrocytes_umap_test_leiden_n15_r0.60.png",
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
        plt.title('Astrocytes Umap Resolution = 2.00 (Leiden)')
        fig.savefig(fname="/home/nolanreardson/Bioinformatics_Env_Windows/Astrocytes_imune_Cells/Clustresults/Astrocytes_umap_test_leiden_n15_r2.00.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

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
        plt.title('Astrocytes Umap Resolution = 1.00 (Leiden)')
        fig.savefig(fname="/home/nolanreardson/Bioinformatics_Env_Windows/Astrocytes_imune_Cells/Clustresults/Astrocytes_umap_test_leiden_n15_r1.00.png",
                    dpi=300,
                    bbox_inches='tight')
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


          

