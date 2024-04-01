# Study clustering across resolutions clustree - LeonorSaude10x
#
# Decide for a good cutoff on resolution
# Requires conversion from .h5ad to .rds 
#
# Daniel Ribeiro, 2023
library("Seurat")
library("R.utils")
library("clustree")

argv <- commandArgs(trailingOnly = FALSE) ##reads the command-line arguments using commandArgs()
base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)) ##sets WD specified in command-line argument
setwd(base_dir)
print(getwd())
options(Ncpus = parallel::detectCores())##sets nº of CPUs 4 parallel processing based on detected nº CPU cores

# Dataset suffix
f_suffix = 'automax'
cell_suffix = '0k'

# Directories
input_dir = "/mnt/c/Users/vasco/Env/Astrocytes_imune_Cells" ##specifies the directory where input data files are located
#input_dir = "Astrocytes_imune_Cells" ##specifies the directory where input data files are located
output_dir = paste0(input_dir, "/Clustresults") ##specifies directory where output files will be saved

# Dataset to analyze
datasets_divided = c('Immune', 'Astrocyte')

# Run clustree
print("Run clustree ...")
clust_prefix = "test_leiden_n15_r"
node_colour = NULL
core_edges = FALSE
layout = "tree"
fig_size = list()   # w x h
fig_size[['Immune']] = c(15, 10)
fig_size[['Astrocyte']] = c(15, 10)

for (d in datasets_divided) ##loop to run clustree for each dataset specified in datasets_divided
{
  filename_in = paste0(input_dir, "/adata_final_", d, "_cca_features.rds")
  print(filename_in)   ##for each dataset, constructs input file path
  adata.seurat = readRDS(filename_in)   ##reads the input data and stores it in the variable adata.seurat
  if (!is.null(node_colour))    ##runs clustree with different options based on node_colour variable
    {
      f_path = paste0(output_dir, "/clustree_plot_", "adata_final_", d, layout, "_", node_colour, ".pdf")
      result = clustree(adata.seurat,
                        layout = layout,
                        prefix = clust_prefix,
                        node_colour = node_colour,
                        use_core_edges = core_edges,
                        return = "plot")

    } else
    {
      f_path = paste0(output_dir, "/clustree_plot_", "adata_final_", d, "_", layout, ".pdf")
      result = clustree(adata.seurat,
                        layout = layout,
                        prefix = clust_prefix,
                        use_core_edges = core_edges,
                        return = "plot")
    } ##The f_path variable holds the output file path, which is constructed based on dataset name and chosen options. 

  ggsave(f_path, plot = result, width = fig_size[[d]][1], height = fig_size[[d]][2], units = "in", limitsize = FALSE) ##Saves the resulting plot using the ggsave() function from the ggplot2 package
}


printf("\n********\n* DONE *\n********\n")
