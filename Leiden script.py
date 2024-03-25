# Cluster tests - LeonorSaude10x
#
# Generate clusters at multiple resolutions.
# To be analyzed in clustree
#
#   Daniel Ribeiro, 2022
import os
import gc
import scanpy as sc

Astrocytes_imune_Cells = "/mnt/c/Users/vasco/Env/Astrocytes_imune_Cells"

### A função já não precisa de n_proc... ###
def start():
    import src.globals
    from src.globals import n_neighbors_final, lineage_resolution_tests 
    from src.globals import data, datasets_divided
    

    for d in ['Astrocyte', 'Immune']:  ##The code iterates over datasets stored in the datasets_divided list
        # Load batch correct data 
        dest = f"{Astrocytes_imune_Cells}/adata_final_{d}_cca_features.h5ad" 
        if os.path.exists(dest): ##checks if the batch-corrected data already exists 
            print("Load batch corrected data...")
            print(dest)
            data[d] = sc.read_h5ad(dest) ##If the file exists, it is loaded using sc.read_h5ad
            keys_to_delete = [key for key in data[d].obs_keys() if key.startswith("test_leiden_n")]
            for key in keys_to_delete:
                del data[d].obs[key] ##alternativa seria data[d].obs.pop(key, None)

        else:
            continue

        ##Setting Up Clustering Parameters
        resolution_tests = lineage_resolution_tests[d]##retrieves the resolution tests specific to a dataset
         
        for i in n_neighbors_final: 
            for r in resolution_tests:
                sc.tl.leiden(adata=data[d], resolution= r, neighbors_key= f'neighbors_{i}',key_added=f'test_leiden_n{i}_r{r:.2f}')
    
            # Free memory
            gc.collect()

        # ## Checkpoint for Maximum Cluster Resolution
        # Here we save a AnnData to assess maximum cluster resolution
        print("After cluster tests:")
        print(data[d])
        # ### Save feature selected version of the data
        print("Save feature selected version of the data...")
        dest = f"{Astrocytes_imune_Cells}/adata_final_{d}_cca_features.h5ad"
        print(dest)
        data[d].write_h5ad(dest, compression='gzip')
        print(data[d])

start()


