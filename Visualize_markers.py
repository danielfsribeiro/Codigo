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
    from globals import H5AD_DIR, Clustresults
    from globals import DATASETS, ENV
    #from helper_functions import random_colors
    
    # STEPS TODO
    # 1. Load literature marker files.
    
    for d in DATASETS:
        #fname = f"{d}5"
        fname = d
        dest = f"{ENV}/{fname}_Markers.csv"
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
            
            #print(adata)
        
        else:
            continue
        
        cell_type_biomarkers = marker_df.groupby('Cell')['Gene'].apply(list).to_dict()

        #print(cell_type_biomarkers)
        #order = marker_df["Cell"].to_list()

        # 3. Visualize DGE markers )
        adata_dge = adata.uns["rank_genes_groups_leiden_fusion"]
        gene_markers = marker_df["Gene"]
        cell_types = marker_df["Cell"].unique()
        vmin = -3
        vmax = 3
        n_genes = 5
        n_groups = adata.obs["leiden_fusion"].cat.categories.size
        x_size = (n_genes * n_groups) * 0.6
        y_size = n_groups

        getrankgenes = sc.get.rank_genes_groups_df(adata=adata,
                                                   group=None,
                                                   key="rank_genes_groups_leiden_fusion",
                                                   pval_cutoff=0.05,
                                                   log2fc_min=0.584,
                                                   log2fc_max=None,
                                                   gene_symbols=None)
        print(getrankgenes)
        getrankgenes.to_csv(f"{ENV}/{d}_rank_genes_filtered.tsv", sep="\t")

        final_genes = {}    # Final list of genes to be plotted in rank_genes_dotplot
        
        # Collect genes in a DataFrame
        for cluster in adata.obs["leiden_fusion"].cat.categories:
            print("Cluster: ", cluster)
            # Filter by log2FC and p-val
            temp_df = sc.get.rank_genes_groups_df(adata=adata,
                                                  group=cluster,
                                                  key="rank_genes_groups_leiden_fusion",
                                                  pval_cutoff=0.05,
                                                  log2fc_min=0.584,     # 1.5x
                                                  log2fc_max=None,
                                                  gene_symbols=None)
            # Filter also by percentage
            mask = temp_df['pct_nz_group'] >= 0.3   # >=30% expression in cluster
            temp_df = temp_df[mask]
            # Sort by log2FC
            temp_df = temp_df.sort_values(by='logfoldchanges', ascending=False)
            print(temp_df)
            print("")

            # Ignore Imm.NA
            if cluster == "Imm.NA":
                continue
            # Collect first 5 genes if there's enough genes
            if temp_df.shape[0] >= 5:
                final_genes[cluster] = temp_df['names'].to_list()[0:5]
            # If not, collect the existing ones, and fill with LogFC high genes the rest of the 5 genes
            else:
                total = 5 - temp_df.shape[0]        # How many genes to fill
                final_genes[f"{cluster}*"] = temp_df['names'].to_list()
                temp_df = sc.get.rank_genes_groups_df(adata=adata,
                                                      group=cluster,
                                                      key="rank_genes_groups_leiden_fusion",
                                                      pval_cutoff=None,
                                                      log2fc_min=0.584,
                                                      log2fc_max=None,
                                                      gene_symbols=None)
                mask = temp_df['pct_nz_group'] >= 0.3   # >=30% expression in cluster
                temp_df = temp_df[mask]
                temp_df = temp_df.sort_values(by='logfoldchanges', ascending=False)
                # Give an indication that the cluster is not OK (*)
                final_genes[f"{cluster}*"].extend(temp_df['names'].to_list()[0:total])
        print(final_genes)

        # How to select some groups only
        mask = adata.obs.loc[:, "leiden_fusion"] != "Imm.NA"
        plotted_adata = adata[mask]
        plotted_adata.obs["leiden_fusion"] = plotted_adata.obs["leiden_fusion"].cat.remove_unused_categories()
        #del plotted_adata.uns["sample_id_colors"]
        #print(plotted_adata)

        # Code for gensbases on Wilcoxon score
        #fig, ax = plt.subplots(figsize=(x_size, y_size))
        fig = sc.pl.rank_genes_groups_dotplot(adata=plotted_adata,
                                              groupby=None,
                                              n_genes=n_genes,
                                              values_to_plot="logfoldchanges",
                                              min_logfoldchange=1,
                                              key="rank_genes_groups_leiden_fusion",
                                              dendrogram=False,
                                              use_raw=None,
                                              show=False,
                                              vmin=-3,
                                              vmax=3,
                                              cmap='bwr',
                                              colorbar_title="log fold change",
                                              return_fig=True)
        
        plt.title(f'{d} Clusters Differential Gene Expression')
        #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        fig.savefig(filename=f"{Clustresults}/{d}_denovo_dotplot.png", dpi=150, bbox_inches='tight')
        plt.close()

        # Code for genes based on DGE>1.5x, adj-pval>0.05, pts>=0.30, sorted by Log2FC
        fig = sc.pl.rank_genes_groups_dotplot(adata=plotted_adata,
                                              groupby=None,
                                              values_to_plot="logfoldchanges",
                                              var_names=final_genes,
                                              key="rank_genes_groups_leiden_fusion",
                                              dendrogram=False,
                                              use_raw=None,
                                              show=False,
                                              vmin=-3,
                                              vmax=3,
                                              cmap='bwr',
                                              colorbar_title="log fold change",
                                              return_fig=True)
        
        plt.title(f'{d} Clusters Differential Gene Expression')
        #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        fig.savefig(filename=f"{Clustresults}/{d}_filtered_denovo_dotplot.png", dpi=150, bbox_inches='tight')
        plt.close()

        # 4. Vizualize literature markers (sc.pl.dotplot) -> dar a lista de genes de marcadores e a lista de tipos de celulas
        n_genes = marker_df.index.size
        x_size = n_genes * 0.4
        
        fig, ax = plt.subplots(figsize=(x_size, 7))
        sc.pl.dotplot(adata=plotted_adata,
                      var_names=cell_type_biomarkers,
                      groupby="leiden_fusion",
                      use_raw=None,
                      categories_order=None,
                      log=False,
                      standard_scale='var',
                      ax=ax,
                      )
        
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        plt.title(f'{d} Cluster Cell Type and Gene Expression')
        fig.savefig(fname=f"{Clustresults}/{fname}_Cell_Type_Gene_expression_dotplot.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

        #sc.pl.dotplot(adata = adata, var_names = marker_df, groupby = cell_types, use_raw=None, log=False, num_categories=0)

        #for cell_type in cell_types:
        #   sc.pl.dotplot(adata, var_names=gene_markers[gene_markers.isin(marker_df.loc[marker_df['Cell'] == #cell_type, 'Gene'])], groupby='leiden_fusion', use_raw=None, log=False, num_categories=0)
        
        # 5. Vizualize literature markers (pl.umap), do a umap for each gene
        #for g in gene_markers:
        #    fig, ax = plt.subplots(figsize=(7, 7))
        #    sc.pl.umap(adata,
        #               use_raw=False,          # gene expression
        #               color=g,
        #               wspace=0.3,
        #               #ncols=3,
        #               cmap="viridis",
        #               neighbors_key='neighbors_15',
        #               show=False,
        #               ax=ax,
        #               vmin = 0,
        #               vmax = 2.0)
        #    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
        #    plt.title(f"{d} '{g}' Expression")
        #    fig.savefig(fname=f"{Clustresults}/{d}_umap_{g}_gene_expression.png", dpi=150, #bbox_inches='tight')
        #    plt.close(fig)

        #cell_types = {}
        #cell_types['M1 macrophage'] = ['CD80', 'CD86']
        #cell_types['M2 macrophage'] = ['MRC1', 'IL-10']


start()
