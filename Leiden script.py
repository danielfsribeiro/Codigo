# Cluster tests - LeonorSaude10x
#
# Generate clusters at multiple resolutions.
# To be analyzed in clustree
#
#   Daniel Ribeiro, 2022
import os
import gc
import scanpy as sc


### A função já não precisa de n_proc... ###
def start():
    from globals import H5AD_DIR
    from globals import lineage_resolution_tests
    
    for d in ['Astrocyte', 'Immune']:   # The code iterates over datasets stored in the datasets_divided list
        # Load batch correct data 
        dest = f"{H5AD_DIR}/adata_final_{d}_cca_features.h5ad" 
        if os.path.exists(dest):    # checks if the batch-corrected data already exists 
            print("Load batch corrected data...")
            print(dest)
            adata = sc.read_h5ad(dest)    # If the file exists, it is loaded using sc.read_h5ad
            # clear old tests
            keys_to_delete = [key for key in adata.obs_keys() if key.startswith("test_leiden_n")]
            for key in keys_to_delete:
                del adata.obs[key]    # alternativa seria data[d].obs.pop(key, None)

        else:
            continue

        # Setting Up Clustering Parameters
        resolution_tests = lineage_resolution_tests[d]      # retrieves the resolution tests specific to a dataset
         
        for r in resolution_tests:
            sc.tl.leiden(adata=adata, resolution=r, neighbors_key='neighbors_15', key_added=f'test_leiden_n15_r{r:.1f}')
            # Free memory
            gc.collect()

        # Checkpoint for Maximum Cluster Resolution
        # Here we save a AnnData to assess maximum cluster resolution
        print("After cluster tests:")
        print(adata)
        # ### Save feature selected version of the data
        print("Save feature selected version of the data...")
        dest = f"{H5AD_DIR}/adata_final_{d}_cca_features.h5ad"
        print(dest)
        adata.write_h5ad(dest, compression='gzip')
        print(adata)


start()
