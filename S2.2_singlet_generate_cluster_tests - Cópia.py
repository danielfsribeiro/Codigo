# Cluster tests - LeonorSaude10x
#
# Generate clusters at multiple resolutions.
# To be analyzed in clustree
#
#   Daniel Ribeiro, 2022
import os
import gc
import scanpy as sc

Astrocytes_imune_Cells = "/mnt/c/Users/vasco/OneDrive/Área de Trabalho/Bioinformatica_IMM/2023-07-16_1809/Astrocytes_imune_Cells"

def start(n_proc=None) -> None:
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
        else:
            continue

        ##Setting Up Clustering Parameters
        resolution_tests = lineage_resolution_tests[d]##retrieves the resolution tests specific to a dataset
        n_tasks = len(n_neighbors_final) * len(resolution_tests) ##???
        print(f"Clustering tasks to execute: {n_tasks}")##calculates the total number of clustering tasks
        if n_proc is None:##checks if n_proc is provided. If not, det. the nº of processes for multiprocessing
            if os.cpu_count() < 4: ##if nº of CPU cores is < than 4, sets n_proc to nº of available CPU cores
                n_proc = os.cpu_count()
            else:
                n_proc = 4 if n_tasks > 4 else n_tasks ##Otherwise, sets n_proc to 4 or nº of tasks, whichever is smaller
        
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


##def leiden_clustering(kwargs: dict):
    ##print(f"Leiden clustering for {kwargs['key_added']} ...")
    ##sc.tl.leiden(**kwargs)
    ##gc.collect()
    ##return kwargs['adata'].obs[kwargs['key_added']].copy()


# main guard required because processes are spawn (compatible with Windows)
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=4)

        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
