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
    from globals import H5AD_DIR, Clustresults
    from globals import DATASETS, ENV
    from helper_functions import random_colors
    

    # STEPS TODO
        # 1. Load literature marker files.
    
    for d in DATASETS:
        
        dest = f"{ENV}/{d}_Markers.csv"
        if os.path.exists(dest):
            print("Load cluster fusion file...")
            print(dest)
            marker_df = pd.read_csv(dest,
                                     sep=';',
                                     header=0,
                                     dtype='category',
                                     index_col=None,
                                     na_filter=False)
            print(marker_df)
        else:
            continue

        # 2. Load entire data (with all genes). File 'raw_norm_annot.h5ad'
        
        dest = f"{H5AD_DIR}/adata_final_{d}_raw_norm_ranked.h5ad"
        if os.path.exists(dest):    # checks if the batch-corrected data already exists
            print("Load batch corrected data...")
            print(dest)
            adata = sc.read_h5ad(dest)
            
            print(adata)
        
        else:
            continue
        
        cell_type_biomarkers = marker_df.groupby('Cell')['Gene'].apply(list).to_dict()

        print(cell_type_biomarkers)

        
        


        # 3. Visualize DGE markers ) 

        adata_dge = adata.uns["rank_genes_groups_leiden_fusion"] 
        gene_markers = marker_df["Gene"]
        cell_types = marker_df["Cell"].unique()


        fig, ax = plt.subplots(figsize=(25, 25))

        sc.pl.rank_genes_groups_dotplot(adata= adata, 
                                        groupby= None, 
                                        values_to_plot= "logfoldchanges", 
                                        min_logfoldchange=1, 
                                        key="rank_genes_groups_leiden_fusion",
                                        use_raw=None, 
                                        show=False,
                                        vmin=-3,
                                        vmax=3,  
                                        cmap='bwr',
                                        ax=ax)
        
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title(f'{d} Clusters Differential Gene Expression')
        fig.savefig(fname=f"{Clustresults}/{d}_denovo_dotplot.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

        # 4. Vizualize literature markers (sc.pl.dotplot) -> dar a lista de genes de marcadores e a lista de tipos de celulas
        
        fig, ax = plt.subplots(figsize=(7, 7))
    

        sc.pl.dotplot(adata = adata,
                      var_names = cell_type_biomarkers, 
                      groupby = "leiden_fusion", 
                      use_raw=None, 
                      log=False,
                      standard_scale ='var', 
                      ax=ax,
                      )  
        
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title(f'{d} Cluster Cell Type and Gene Expression')
        fig.savefig(fname=f"{Clustresults}/{d}_Cell_Type_Gene_expression_dotplot.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

        #sc.pl.dotplot(adata = adata, var_names = marker_df, groupby = cell_types, use_raw=None, log=False, num_categories=0)

        #for cell_type in cell_types:
            #sc.pl.dotplot(adata, var_names=gene_markers[gene_markers.isin(marker_df.loc[marker_df['Cell'] == #cell_type, 'Gene'])], groupby='leiden_fusion', use_raw=None, log=False, num_categories=0)
        




        # 5. Vizualize literature markers (pl.umap), do a umap for each gene
        

        fig, ax = plt.subplots(figsize=(7, 7))
        colors = random_colors(len(adata_dge.cat.categories.to_list()))

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
        plt.title(f'{d} Differential Gene Expression')
        fig.savefig(fname=f"{Clustresults}/{d}_umap_differential_gene_expression.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

        #cell_types = {}
        #cell_types['M1 macrophage'] = ['CD80', 'CD86']
        #cell_types['M2 macrophage'] = ['MRC1', 'IL-10']

       
        
        

       


start()
