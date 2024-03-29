import os
import gc
import scanpy as sc
import pandas as pd
import numpy as np

### Defines
Good_clusters_r2 = ['1', '10', '15', '24', '30']      # .obs['test_leiden_n15_r2']
Good_clusters_r1 = ["6", "7", "9"]
Bad_clusters_r1 = ["0", "1", "2", "3", "4,", "5", "8"]


# Convert old cluter nr into new cluster nr
dict_good_cells = {'1': '2.1',
                   '10': '2.10',
                   "15": "2.15",
                   "24": "2.24",
                   "30": "2.30",
                   }
dict_bad_cells = {'0': '1.0',
                  '1': '1.1',
                  '2': '1.2',
                  '3': '1.3',
                  '4': '1.4',
                  '5': '1.5',
                  '8': '1.8',
                  }

### Load data
# Use adata_cca_features.h5ad --> must contain info about r1 and r2

Astrocytes_imune_Cells = "/mnt/c/Users/vasco/Env/Astrocytes_imune_Cells"


def start():
    from globals import data, data2

    for d in ['Immune']:  ##The code iterates over datasets stored in the datasets_divided list
        # Load batch correct data 
        dest = f"{Astrocytes_imune_Cells}/adata_final_{d}_cca_features.h5ad" 
        if os.path.exists(dest):    ##checks if the batch-corrected data already exists 
            print("Load batch corrected data...")
            print(dest)
            data[d] = sc.read_h5ad(dest)
            adata = data[d]     ##If the file exists, it is loaded using sc.read_h5ad             
        else:
            continue
          
        mask_r1 = adata.obs["test_leiden_n15_r1.00"].isin(Good_clusters_r1)
        mask_r2 = adata.obs["test_leiden_n15_r2.00"].isin(Good_clusters_r2)
        Goodcells = adata.obs[mask_r1 & mask_r2]
        Badcells = adata.obs[~(mask_r1 & mask_r2)]
        cell_names_good = Goodcells.index.tolist()
        cell_names_bad = Badcells.index.tolist()

        for d in ['Immune']:   
            dest2 = f"{Astrocytes_imune_Cells}/adata_final_{d}_cca_features.h5ad" 
            if os.path.exists(dest2): ##checks if the batch-corrected data already exists
                print("Load batch corrected data...")
                print(dest2)
                data2[d] = sc.read_h5ad(dest2)
                adata2 = data2[d] ##If the file exists, it is loaded using sc.read_h5ad
                adata2.obs['Fusion_cluster'] = None
            else:
                continue

            adata2 = adata2[cell_names_good, :].obs['Fusion_cluster'] = Goodcells
          
            ## Convert cluster names
            for key, value in dict_good_cells.items():
                adata2[adata2.obs['Fusion_cluster'].isin([key])].obs['Fusion_cluster'] = value
            
            print(adata2.obs['Fusion_cluster'])
  
          
start()
