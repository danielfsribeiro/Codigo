# Calculate cell fractions in clusters - LeonorSaude10x
# 1. Calculate if some cluster is enriched for a particular sample
# 2. Calculate if some cluster is enriched for a particluar cell type seen by scMCA.
# 3. Plot fractions per cluster
#
#
# Daniel Ribeiro, 2023
import os
import gc
import re
import numpy as np
import pandas as pd

import scanpy as sc
import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from statsmodels.stats.multitest import multipletests ##used for multiple hypothesis testing correction


def sample_factions(adata: sc.AnnData,
                    dataset_type: str,
                    key: str,
                    condition: str = 'injury_condition') -> pd.DataFrame: 
    ##calculates the fractions of cells belonging to different clusters for each sample (condition) in the scRNA-seq datase

    """
    adata: `sc.AnnData`
        Dataset containing cluster information
    key: `str`
        .obs key to look at
    condition: `str`
        String to choose from injury_day or injury_condition list of samples
    """
    from globals import injury_condition, injury_day, injury_region, collection_region, \
        injury_region_no_central, collection_region_no_central

    samples = [] ##creates an empty list samples to store the selected samples
    prog = re.compile(r'^injured', flags=re.IGNORECASE)
    if condition == 'injury_condition': ##condition parameter is checked
        samples = list(injury_condition) ##corresponding samples selected based on dataset type (dataset_type) 
        if dataset_type.find('uinj') != -1:##checks if dataset_type contains the string 'uinj'
            samples = [s for s in samples if prog.search(s) is None]##If+,filt spl names in spl lst=RE.pattern
            ##The same filtering logic is applied for other condition values
    elif condition == 'injury_day':
        samples = list(injury_day)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'injury_region':
        samples = list(injury_region)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'collection_region':
        samples = list(collection_region)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'injury_region_no_central':
        condition = 'injury_region'
        samples = list(injury_region_no_central)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'collection_region_no_central':
        condition = 'collection_region'
        samples = list(collection_region_no_central)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    else:
        raise ValueError("'condition' must be either a string 'injury_condition', 'injury_day', 'injury_region', 'collection_region'.") ##raises a ValueError if an invalid condition value is provided
    
    ##Calculates the cell fractions for each cluster and sample
    # Subset adata 
    ##selects rows (cells) where sample condition matches 1 of selected samples
    adata2 = adata[adata.obs[condition].isin(samples), :].copy()
    
    ##creates new AnnData object by subsetting the original adata object based on the selected samples defined in the samples list. Only keeps the cells that match specified condition (sample name).isin(samples) function filters cells where the condition matches any of sample names in the samples list. Then, .copy() creates a copy of this subset to avoid modifying the original adata object
    # Initiate matrix
    ##initiates a matrix table_frac filled with zeros rows represent clusters, and columns represent samples
    n_cols = len(samples) + 1 ## nº of selected samples +1, will be additional column "total in cluster"
    n_rows = len(adata2.obs[key].cat.categories) + 2## nº clusters in adata2.obs specified key DataFrame (cell annotations) +2, for "frac sample in total" and "total in sample" rows
    table_frac = np.zeros((n_rows, n_cols))##empty 2D np array, dim(n_rows,n_cols) to store calc. fractions
    total_cluster = adata2.obs[key].value_counts().loc[adata2.obs[key].cat.categories].values.tolist() 
    ##counts the number of cells in each cluster and stores them in total_cluster
   
    # Reorder samples logically
    sample_counts = adata2.obs[condition].value_counts().loc[samples].copy()
    ##counts nº of cells in each selected sample based on condition col of adata2.obs DataFrame
    total_sample = sample_counts.values.tolist()##counts are converted to a list and stored in total_sample
    
    # Calculate fractions
    ##Nested loop used to calculate the fraction of cells in each cluster (row) belonging to each sample (col).
    for i in range(len(total_cluster)):
        for j in range(len(sample_counts.index)):
            # Cells of cluster i that belong to sample j
            n_cells = adata2[(adata2.obs[key] == adata2.obs[key].cat.categories[i]) & \
                             (adata2.obs[condition] == sample_counts.index[j])].n_obs
            table_frac[i, j] = n_cells / total_cluster[i]
   
    # Calculate sample fractions in total and fill totals
    ##Calc. fraction of cells in each spl in entire dataset (in all clusters), fils totals 4 each selected spl
    row_i = len(total_cluster)# fraction in total##set to the index where the "frac sample in total" row should be placed in the table_frac matrix
    for j in range(len(sample_counts.index)):##iterate each selected spl,4 each, calc fract cells out Tnºcells
        table_frac[row_i, j] = adata2[adata2.obs[condition] == sample_counts.index[j]].n_obs / adata2.n_obs
        table_frac[row_i + 1, j] = total_sample[j]##result stored in table_frac matrix at corresponding posit.
        ##fills in total nº of cells for each selected spl in the last row of the table_frac matrix
    
    # Fill cluster totals
    ##Fills in the total number of cells for each cluster in the last column of the table_frac matrix
    col_j = len(total_sample)#totals per cluster##col_j set to index where"total in cluster"col should be placed in the table_frac matrix.
    for i in range(len(total_cluster)):##iterat each clust.i 4 each fills nºTcell in last col table_frac matrix
        table_frac[i, col_j] = total_cluster[i]
    # Total cells in dataset##fill in Tnº cells in entire dataset in bottom-right cell of table_frac matrix
    table_frac[n_rows - 1, n_cols - 1] = adata2.n_obs
    table_frac[n_rows - 2, n_cols - 1] = adata2.n_obs / adata2.n_obs

    # Export table
    ##Sets appropriate row and column names for the table_frac matrix and converts it into a pandas DataFrame
    row_names = adata2.obs[key].cat.categories.to_list()
    row_names.append('frac sample in total')
    row_names.append('total in sample')
    col_names = samples
    col_names.append('total in cluster')
    table_frac = pd.DataFrame(table_frac, index=row_names, columns=col_names)

    return table_frac


def sample_counts(adata: sc.AnnData,##calc & return DF(table_counts) representing cell cnt 4 each cluster & spl
                  dataset_type: str,##str indicating type of dataset (e.g.'uinj') used to determine list of spl
                  key: str, ##represents the .uns key containing sample fractions in clusters x samples format
                  condition: str = 'injury_condition') -> pd.DataFrame:##str that allows choosing from a set of predefined sample lists, such as 'injury_condition', 'injury_day', 'injury_region', 'collection_region', etc. The default value is 'injury_condition'.
    """
    adata: `sc.AnnData`
        Dataset containing cluster information
    key: `str`
        .uns key containing sample fractions (clusters x samples)
    condition: `str`
        String to choose from injury_day or injury_condition list of samples
    """
    from globals import injury_condition, injury_day, injury_region, collection_region, \
        injury_region_no_central, collection_region_no_central##imports vars from src.globals

    samples = []
    prog = re.compile(r'^injured', flags=re.IGNORECASE)##reg expres. patt. prog created to match spl names start with'injured'
    if condition == 'injury_condition':##sets the samples list based on the provided condition
        samples = list(injury_condition)##If condition is 'injury_condition',sets spl to list injury_condition
        if dataset_type.find('uinj') != -1:##checks if dataset_type contains the string 'uinj'
            samples = [s for s in samples if prog.search(s) is None]##If+,filt spl names in spl lst=RE.pattern
            ##The same filtering logic is applied for other condition values 
    elif condition == 'injury_day':
        samples = list(injury_day)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'injury_region':
        samples = list(injury_region)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'collection_region':
        samples = list(collection_region)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'injury_region_no_central':##add condition modifies spl lst based on prvided dataset_type
        condition = 'injury_region'
        samples = list(injury_region_no_central)##update spl list with the corresponding predefined lists 
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'collection_region_no_central':##add condition modifie spl lst based prvided dataset_type
        condition = 'collection_region'
        samples = list(collection_region_no_central)##update spl list with the corresponding predefined lists
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    else:
        raise ValueError("'condition' must be either a string 'injury_condition', 'injury_day', 'injury_region', 'collection_region'.")##if condition do not match predef conditions valuerror raised

    # Subset adata
    ##Creates a new AnnData object called adata2 by subsetting the original adata object based on the selected samples defined in the samples list. It only keeps the cells (rows) in adata that match the specified condition (sample name) in the condition variable. The .isin(samples) function filters the rows where the condition matches any of the sample names in the samples list. Then, .copy() creates a copy of this subset to avoid modifying the original adata object.
    adata2 = adata[adata.obs[condition].isin(samples), :].copy()
    # Initiate matrix
    n_cols = len(samples) + 1##nº of selected spls+1,there will be an add col 4 "total in cluster"
    n_rows = len(adata2.obs[key].cat.categories) + 1##nº clusters in specified key of adata2.obs DataFrame +1

    table_counts = np.zeros((n_rows, n_cols))##empty np array called table_counts created dim(n_rows, n_cols)
    ##This array will be used to store the calculated cell counts 
    total_cluster = adata2.obs[key].value_counts().loc[adata2.obs[key].cat.categories].values.tolist()   # Order clusters
    ##The variable total_cluster is calculated to store the count of cells in each cluster present in adata2. It first counts the number of cells belonging to each cluster (using value_counts()), considering only the categories (clusters) present in the key column of the adata2.obs DataFrame. The counts are then converted to a list.

    # Reorder samples logically
    sample_counts = adata2.obs[condition].value_counts().loc[samples].copy()##stores the count of cells for each selected sample in adata2. It counts the number of cells belonging to each sample based on the condition column of the adata2.obs DataFrame
    total_sample = sample_counts.values.tolist()##counts are converted to a list and stored in total_sample.
    
    # Calculate fractions
    for i in range(len(total_cluster)):##Nested loop used to calc cell cnts 4 each cluster (row) & spl (column)
        for j in range(len(sample_counts.index)):##4 each cluster i & spl j, calc nº cells from clust i & spl j
            # Cells of cluster i that belong to sample j
            n_cells = adata2[(adata2.obs[key] == adata2.obs[key].cat.categories[i]) & \
                             (adata2.obs[condition] == sample_counts.index[j])].n_obs
            table_counts[i, j] = n_cells##result stored in the table_counts matrix

    # Total sample counts
    ##calculates the T nº cells for each selected spl and fills them in the last row of the table_counts matrix
    row_i = len(total_cluster)#fraction in T##row_i set->index whr"total in sample"row shoud be table_counts
    for j in range(len(sample_counts.index)):##loop itrate each sel spl j & 4 each fils Tnºcel of that spl
        table_counts[row_i, j] = total_sample[j]##result stored in table_counts mtx at corresponding position
    
    # Fill cluster totals
    ##Fills in the total number of cells for each cluster in the last column of the table_counts matrix
    col_j = len(total_sample)#totals per cluster##set to idx whr"total in cluster"col shoud be in table_counts
    for i in range(len(total_cluster)):##loop iterates each cluster i, & 4 each, fils Tnºcell in that cluster 
        table_counts[i, col_j] = total_cluster[i]##result stored in table_counts mtx at corresponding position

    # Total cells in dataset
    ##last cell of table_counts matrix (bottom-right cell) filled with Tnºcell in entire dataset (adata2.n_obs)
    table_counts[n_rows - 1, n_cols - 1] = adata2.n_obs

    # Export table
    #Sets appropriate row and column names for the table_counts matrix and converts it into a pandas DataFrame
    row_names = adata2.obs[key].cat.categories.to_list()
    row_names.append('total in sample')##clusters in key col of adata2.obs DF & add row 4 "total in sample"
    col_names = samples##contains the selected samples along with an additional column name "total in cluster"
    col_names.append('total in cluster')
    table_counts = pd.DataFrame(table_counts, index=row_names, columns=col_names)

    return table_counts


def plot_stacked_plots(data: pd.DataFrame,
                       dataset_name: str,
                       alternate_sample_names=None
                       ) -> None:
    from globals import fractions_dir
    from globals import proportion_colors
    
    clusters = [i for i in range(len(data.index))]##list created, containing int from 0 to len of index DF This list is used as x-axis ticks in the subsequent plot
    colors = [proportion_colors[c] for c in data.columns]##list created by mapping the colors from the proportion_colors dictionary to the column names of the input DataFrame data
    bottom = np.zeros(len(clusters))##array initialized with zeros and has the same length as the clusters list. It will be used to keep track of the bottom positions of each bar in the stacked plot

    with PdfPages(fractions_dir + f"/{dataset_name}_stacked_plot.pdf", keep_empty=True) as pdf:
        ##PDF file is created using the PdfPages class from the matplotlib.backends.backend_pdf module, file saved in the directory specified by fractions_dir with a filename based on the dataset_name
        fig, ax = plt.subplots(figsize=(len(data.index) * 0.5, 7))##new fig & axes created using plt.subplots(). Figsize parameter sets the size of the figure based on the number of rows in the DataFrame data
        ax.grid(False)## grid lines are turned off on the plot.
        ##loop iterates through the columns of the DataFrame data. For each column, a bar plot is created using the bar method of the ax object. The x-coordinates of the bars are specified by clusters, the heights are taken from the corresponding column of data, and the bottom parameter is used to stack the bars on top of each other
        for i, b in enumerate(data.columns):
            ax.bar(clusters,
                   height=data.iloc[:, i].to_numpy().flatten(),
                   bottom=bottom,
                   color=colors[i],
                   width=0.6,
                   yerr=0,
                   label=b)
            # Next iteration bars will be plotted on top of the previous one
            bottom += data.iloc[:, i].to_numpy().flatten()

        ax.set_xticks(clusters)##The x-axis ticks are set using ax.set_xticks(clusters)
        ##If alternate_sample_names is provided (a list of alternate names for the x-axis ticks), the x-axis labels are set to those names with a 45-degree rotation and alignment to the right.
        if alternate_sample_names is not None:
            ax.set_xticklabels(alternate_sample_names, rotation=45, ha='right', rotation_mode='anchor')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        pdf.savefig(dpi=300, bbox_inches='tight')##plot saved to PDF using pdf.savefig(). dpi sets resolution, and bbox_inches='tight' ensures that the saved figure does not have any unnecessary white spaces. 
        plt.close(fig)##closes the figure to free up memory

##Chisquare function from the scipy.stats module to perform a chi-square test. Takes a dictionary kwargs as input, which contains the necessary arguments for the chisquare function. Chi-square test is used to compare the observed frequencies in a cluster to the expected frequencies in the cluster given the frequencies in the entire dataset. The data for the chi-square test is in the form of a contingency table, where rows represent clusters and columns represent conditions (e.g., injury_condition, injury_day, etc.). The function returns the results of the chi-square test.

def chisquare_mp(kwargs: dict):
    # Given the cluster size, calculate the expected cell dist if they respected dataset dist
    # Chi square test
    # m (rows) x n (cols) -> observed in cluster vs expected in cluster given the frequencies in dataset
    #          | cond1 | cond2 | cond3 | ...
    # -------------------------------------
    # cluster  | obs1  | obs2  | obs3  |
    # ------------------------------------
    # expected | exp1  | exp2  | exp3  |
    from scipy.stats import chisquare
    return chisquare(**kwargs)

##This function uses the fisher_exact function from the FisherExact module to perform a Fisher's exact test. It takes a dictionary kwargs as input, which contains the necessary arguments for the fisher_exact function. The Fisher's exact test is used to compare the observed frequencies in a cluster to the observed frequencies in the entire dataset. The data for the Fisher's exact test is also in the form of a contingency table, similar to the chi-square test. The function returns the results of the Fisher's exact test.

def fishertest_mp(kwargs: dict):
    # m (rows) x n (cols) -> observed in cluster vs observed in dataset
    #         | cond1 | cond2 | cond3 | ...
    # -------------------------------------
    # cluster | a     | b     | c     |
    # ------------------------------------
    # dataset | d     | e     | f     |
    from FisherExact import fisher_exact
    return fisher_exact(**kwargs)

##The start function is the main function that initiates the processing of data and statistical tests. It loads necessary modules and data from the src.globals module. The function takes an optional argument n_proc, which determines the number of processes for multiprocessing

def start(n_proc=None) -> None:
    import globals
    from globals import checkpoint_dir, fractions_dir
    from globals import data, datasets_divided, lineage_resolution_final, n_neighbors_final
    
    # Load ranked data
    #for d in list(datasets_divided.keys()):
    #    datasets_divided[f'{d}_uinj'] = datasets_divided[d]
    for d in datasets_divided:##loop iterates datasets_divided dict & loads gene rank data from H5AD files
        dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked.h5ad"
        if os.path.exists(dest):
            print("Load gene rank data...")
            print(dest)
            data[d] = sc.read_h5ad(dest)##The loaded data is stored in the data dictionary
        else:
            continue
        
        # Calculate fractions and totals
        ##Inside another loop, the function calculates fractions and totals for each condition based on the clustering results. It uses functions sample_factions and sample_counts, passing the data, dataset_type, key (clusters), and condition as arguments.The calculated fractions and counts are saved as CSV files and stored in the AnnData object.
        for res in lineage_resolution_final[d]:
            clusters = f'leiden_n{n_neighbors_final[0]}_r{res}'
            for condition in ['injury_condition', 'injury_day', 'injury_region', 'collection_region', \
                              'injury_region_no_central', 'collection_region_no_central']:
                print(f"Calculate statistics for {condition}...")
                # Save fraction
                table_frac = sample_factions(data[d], dataset_type=d, key=clusters, condition=condition)
                dest = f"{fractions_dir}/{d}_final_sample_frac_{clusters}_{condition}.txt"
                print(f"Output: {dest}")
                table_frac.to_csv(dest, sep='\t')
                # Update AnnData
                data[d].uns[f'sample_frac_{clusters}_{condition}'] = table_frac
                
                # Save totals
                table_count = sample_counts(data[d], dataset_type=d, key=clusters, condition=condition)
                dest = f"{fractions_dir}/{d}_final_sample_counts_{clusters}_{condition}.txt"
                print(f"Output: {dest}")
                table_count.to_csv(dest, sep='\t')
                # Update AnnData
                data[d].uns[f'sample_totals_{clusters}_{condition}'] = table_count
                
                # Do chi-square and fisher tests
                ##The data for these tests is taken from the table_count DataFrame, which contains observed counts for each cluster and condition. chisquare_mp and fishertest_mp functions are used for the tests, and the results are saved in a DataFrame df, which is then stored as a CSV file.The code calculates the p-values and adjusted p-values (FDR correction) for both chi-square and Fisher's exact tests. The results are stored in the df DataFrame and saved as a CSV file.
                
                expected_freq = table_frac.iloc[-2, :-1].values  # Expected freqeuncies observed in data set
                chi_test_result = []
                fisher_result = []
                seed = 2015854237
                
                # Multiprocessing
                from multiprocessing.pool import Pool
                n_tasks = len(table_count.index[:-1])
                print(f"Tasks to execute: {n_tasks}")
                if n_proc is None:
                    if os.cpu_count() < 4:
                        n_proc = os.cpu_count()
                    else:
                        n_proc = 4 if n_tasks > 4 else n_tasks
                with Pool(processes=n_proc) as p:
                    # Chi square test
                    print("Chi square test...")
                    args = [{'f_obs': table_count.iloc[row, :-1].values,
                             'f_exp': [freq * table_count.iloc[row, -1] for freq in expected_freq]
                             }
                            for row in range(len(table_count.index[:-1]))]
                    #print(args)
                    # Returns: reject, pvals_corrected, alphacSidak, alphacBonf
                    chi_test_result = list(p.map(chisquare_mp, args))
                    gc.collect()
                    padj_chi = multipletests([chi[1] for chi in chi_test_result], method='fdr_bh')
                    padj_chi = padj_chi[1].tolist()
                    
                    # Fisher test
                    print("Fisher test...")
                    args = [{'table': [table_count.iloc[row, :-1].values.tolist(), table_count.iloc[-1, :-1].values.tolist()],
                             'seed': seed,
                             'simulate_pval': True,
                             'replicate': 100000,
                             'workspace': 100000000
                             }
                            for row in range(len(table_count.index[:-1]))]
                    # Returns pvalue
                    fisher_result = list(p.map(fishertest_mp, args))
                    gc.collect()
                    padj_fisher = multipletests([f for f in fisher_result], method='fdr_bh')
                    padj_fisher = padj_fisher[1].tolist()
                    
                    df = {'chi.pval': [chi[1] for chi in chi_test_result],
                          'chi.adjpval': padj_chi,
                          'fisher.pval': fisher_result,
                          'fisher.adjpval': padj_fisher}
                    df = pd.DataFrame(df, index=table_count.index[:-1])
                    dest = f"{fractions_dir}/{d}_final_proportion_statistics_{condition}.txt"
                    print(f"Output: {dest}")
                    df.to_csv(dest, sep='\t')
                gc.collect()
                

                # Plot stacked barplots
                ##The function generates stacked barplots based on the plot_data DataFrame. It calls the plot_stacked_plots function, passing plot_data, a dataset name, and alternate sample names as arguments. The generated plots are saved as a PDF file
                print("\nPlot stacked barplots...")
                plot_data = data[d].uns[f'sample_frac_{clusters}_{condition}'].iloc[:-2, :-1].copy()
                plot_data.loc['Expected', :] = expected_freq
                #alternate_names = data[d].uns[f'annot_{clusters}'].loc[:, 'cell_lineage'].values.tolist()
                #alternate_names = [f'{i}: {n}' for i, n in enumerate(alternate_names)]
                alternate_names = data[d].uns[f'sample_frac_{clusters}_{condition}'].index[:-2].to_list()
                alternate_names.append('Expected')
                plot_stacked_plots(plot_data,
                                   dataset_name=f'{d}_final_{clusters}_{condition}',
                                   alternate_sample_names=alternate_names)
                # With dendrogram order
                plot_data = data[d].uns[f'sample_frac_{clusters}_{condition}'].iloc[:-2, :-1].copy()
                plot_data = plot_data.loc[data[d].uns[f'dendrogram_{clusters}']['dendrogram_info']['ivl'].tolist(), :]
                plot_data.loc['Expected', :] = expected_freq
                alternate_names = data[d].uns[f'dendrogram_{clusters}']['dendrogram_info']['ivl'].tolist()
                #alternate_names = [f'{plot_data.index[i]}: {alternate_names.iloc[i]}' for i in range(len(plot_data.index))]
                alternate_names.append('Expected')
                plot_stacked_plots(plot_data,
                                   dataset_name=f'{d}_final_{clusters}_{condition}_dendrogram',
                                   alternate_sample_names=alternate_names)
        
        # Save AnnData
        ##After all processing is completed, the function saves the updated data dictionary (AnnData objects) back to H5AD files using the write_h5ad() method from the scanpy library.
        print("Save AnnData...")
        dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked.h5ad"
        data[d].write_h5ad(dest, compression='gzip')


# main guard required because processes are spawn (compatible with Windows)
##The code block ensures that multiprocessing is used by setting the start method for Windows compatibility. It sets the number of processes for multiprocessing to n_proc=8,
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=8)

        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
