# Gene Set enrichment analysis (GSEA) - LeonorSaude10x
#
# Test pre-ranked GSEA with multiple types of rankings:
# 1. Score from Wilcoxon rank-sum test
# 2. Gene enrichment score calculated in S4.1
# 3. t-stat calculated by limma in S6.3
#
# Daniel Ribeiro, 2023
import os
import gc

import pandas as pd
import scanpy as sc




def prerank_gsea(kwargs: dict) -> None:
    import gseapy as gp
    from globals import gene_sets
    
    seed = 2015854237
    d = kwargs['dataset']
    c = kwargs['cluster']
    
    print(f"Calculate for {d} {c}")
    kwargs['scores_df'].sort_values(by='scores', ascending=False, inplace=True)
    kwargs['scores_df'].set_index('names', verify_integrity=True, inplace=True)
                
    # Build prerank df
    # Substitute mouse gene names by human gene names
    prerank_df = pd.DataFrame()
    prerank_df['gene'] = kwargs['converted_genes_df'].loc[kwargs['scores_df'].index[kwargs['scores_df'].index.isin(kwargs['converted_genes_df'].index)]]
    prerank_df.drop_duplicates(subset='gene', inplace=True)
    prerank_df['score'] = kwargs['scores_df'].loc[prerank_df.index, 'scores']
    prerank_df.set_index('gene', verify_integrity=True, inplace=True)
    
    # Prerank GSEA
    # GSEApy supports pandas DataFrame, Series, with gene names as index.
    # Supports text file without any header, just 'gene_name' 'score' per row
    # Run GSEA per gene set
    for k, v in gene_sets.items():
        for gs in v:
            print(f"Analyzing {d} {c} {k} - {gs}...")
            dest = kwargs['gsea_dir'] + k.replace(':', '_') + f"/{d}_{c}_{gs}"
            pre_res = gp.prerank(rnk=prerank_df,
                                 gene_sets=gs,
                                 threads=kwargs['threads'],
                                 min_size=10,   # 5
                                 max_size=2000,  # [500-1000]
                                 permutation_num=1000,  # reduce number to speed up testing
                                 graph_num=50,   # 20
                                 outdir=dest,
                                 seed=seed,
                                 verbose=True)
            #print(pre_res)
        
        # Free memory
        gc.collect()


def gsea_wilcoxon(dataset: str,
                  selected_genes: pd.DataFrame,
                  h5adadata: sc.AnnData,
                  groupby: str,     # key in .obs
                  unskey : str,     # key in .uns
                  gsea_dir: str,
                  n_proc: int) -> None:
    
    clusters = h5adadata.obs[groupby].cat.categories.to_list()
    for c in clusters:
        getrankgenes = sc.get.rank_genes_groups_df(adata=h5adadata,
                                                    group=c,
                                                    key= unskey)
        print(c, ':')
        print(getrankgenes)
        print()

        # Run GSEA for each cluster
        print("Run pre-ranked GSEA...")
        if n_proc is None:
            n_proc = os.cpu_count() - 1
        args = {'scores_df': getrankgenes.loc[:, ["names", "scores"]].copy(),
                'converted_genes_df': selected_genes,
                'dataset': dataset,
                'cluster': c,
                'gsea_dir': f'{gsea_dir}/wilcoxon_',
                'threads': n_proc,
                }
        prerank_gsea(args)
        
    # Free memory
    gc.collect()


def start(n_proc=None) -> None:
    import mygene
    import time
    import datetime
   
    from globals import H5AD_DIR, GSEA_DIR
    from globals import DATASETS
    
    
    # Load socres
    t1 = time.time()
    for d in DATASETS:
        dest = f"{H5AD_DIR}/adata_final_{d}_raw_norm_ranked.h5ad"
        if os.path.exists(dest):    # checks if the batch-corrected data already exists
            print("Load feature data...")
            print(dest)
            adata = sc.read_h5ad(dest)
            
            print(adata)
            
        else:
            continue
        
        # Convert mouse genes to humam
        t2 = time.time()
        mg = mygene.MyGeneInfo()
        converted = mg.querymany(adata.var_names, scopes='symbol,alias', species='human', fields='ensembl.gene,symbol', as_dataframe=True)
        converted.dropna(axis=0, subset='symbol', inplace=True)
        converted = converted['symbol']
        
        # Marker genes
        print("Calculate GSEA for ranked genes (wilcoxon scores)...")
        gsea_wilcoxon(dataset=d,
                      selected_genes=converted,
                      unskey= "rank_genes_groups_leiden_fusion",
                      h5adadata = adata,
                      groupby= "leiden_fusion",
                      gsea_dir=GSEA_DIR,
                      n_proc=n_proc)
        
        print(f"Time for {d}: {datetime.timedelta(seconds=(time.time()-t2))}")
        
    print(f"Time total: {datetime.timedelta(seconds=(time.time()-t1))}")


# main guard required because processes are spawn (compatible with Windows)
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=None)

        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
