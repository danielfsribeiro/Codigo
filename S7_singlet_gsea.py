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
import matplotlib
matplotlib.use('cairo')


def prerank_gsea(kwargs: dict) -> None:
    import gseapy as gp
    from globals import gene_sets
    
    seed = 2015854237
    d = kwargs['dataset']

    neigh = kwargs['neighbor']
    res = kwargs['resolution']
    c = kwargs['cluster']
    
    print(f"Calculate for {d} leiden_n{neigh}_r{res}_c{c}")
    kwargs['scores_df'].sort_values(by=f'{c}.scores', ascending=False, inplace=True)
    kwargs['scores_df'].set_index(f'{c}.names', verify_integrity=True, inplace=True)
                
    # Build prerank df
    prerank_df = pd.DataFrame()
    prerank_df['gene'] = kwargs['converted_genes_df'].loc[kwargs['scores_df'].index[kwargs['scores_df'].index.isin(kwargs['converted_genes_df'].index)]]
    prerank_df.drop_duplicates(subset='gene', inplace=True)
    prerank_df['score'] = kwargs['scores_df'].loc[prerank_df.index, f'{c}.scores']
    prerank_df.set_index('gene', verify_integrity=True, inplace=True)
    
    # Prerank GSEA
    # GSEApy supports pandas DataFrame, Series, with gene names as index.
    # Supports text file without any header, just 'gene_name' 'score' per row
    # Run GSEA per gene set
    for k, v in gene_sets.items():
        for gs in v:
            print(f"Analyzing {d} leiden_n{neigh}_r{res}_c{c} {k} - {gs}...")
            dest = kwargs['gsea_dir'] + k.replace(':', '_') + f"/{d}_leiden_n{neigh}_r{res}_c{c}_{gs}"
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
                  neighbor: int,
                  resolution: float,
                  n_hvg: int,
                  ranked_genes_dir: str,
                  gsea_dir: str,
                  n_proc: int) -> None:
    dest = f"{ranked_genes_dir}/{dataset}_final_marker_genes_leiden_n{neighbor}_r{resolution}_{n_hvg}.txt"
    print("Load marker gene data...")
    print(dest)
    df_marker = pd.read_csv(dest, sep='\t', index_col=0, header=0)
    
    # Find clusters
    clusters = list(set([s.split('.')[0] for s in df_marker.columns.to_list()]))
    clusters.sort(key=lambda x: int(x))
    
    # Run GSEA for each cluster
    print("Run pre-ranked GSEA...")
    n_tasks = len(clusters)
    print(f"Tasks to execute: {n_tasks}")
    if n_proc is None:
        n_proc = 4
    for c in clusters:
        args = {'scores_df': df_marker.loc[:, [f'{c}.names',
                                               f'{c}.scores',
                                               f'{c}.logfoldchanges',
                                               f'{c}.pvals',
                                               f'{c}.pvals_adj',
                                               f'{c}.pts']].copy(),
                'converted_genes_df': selected_genes,
                'dataset': dataset,
                'neighbor': neighbor,
                'resolution': resolution,
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
    import src.globals   # noqa:F401
    from src.globals import checkpoint_dir, annotation_dir, marker_genes_dir, gsea_dir
    from src.globals import data, datasets_divided, lineage_resolution_final, n_neighbors_final, n_hvg_spi
    
    neigh = n_neighbors_final[0]
    # Load socres
    t1 = time.time()
    for d in datasets_divided:
        dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked.h5ad"
        if os.path.exists(dest):
            print("Load gene rank data...")
            print(dest)
            data[d] = sc.read_h5ad(dest)
        else:
            continue
        
        # Convert mouse genes to humam
        t2 = time.time()
        mg = mygene.MyGeneInfo()
        converted = mg.querymany(data[d].var_names, scopes='symbol,alias', species='human', fields='ensembl.gene,symbol', as_dataframe=True)
        converted.dropna(axis=0, subset='symbol', inplace=True)
        converted = converted['symbol']
        
        # Marker genes
        print("Calculate GSEA for ranked genes (wilcoxon scores)...")
        gsea_wilcoxon(dataset=d,
                      selected_genes=converted,
                      neighbor=neigh,
                      resolution=lineage_resolution_final[d][0],
                      n_hvg=n_hvg_spi,
                      ranked_genes_dir=marker_genes_dir,
                      gsea_dir=gsea_dir,
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
