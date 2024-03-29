#!/bin/bash
conda activate py39

# Exececute a file at a time
# TODO: What does "Leiden script.py" do? Maybe give it a more meaningful name
# Creates cluslters at multiple resolutions
FILE="Leiden script.py"
python "$FILE"
if [[ $? -eq 0 ]]
then
  echo "Success '$FILE'"
  echo "Script ran successfuly."
else
  echo "Error '$FILE'"
  echo "Script did not complete. Exit all analysis."
  return 1
fi

# TODO
# Vizualize clusters. Didn't find code for Immune cells
# TODO: Harmonize name and make it general for all populations
FILE="Astrocytes_GoodResolution.py"
python "$FILE"
if [[ $? -eq 0 ]]
then
  echo "Success '$FILE'"
  echo "Script ran successfuly."
else
  echo "Error '$FILE'"
  echo "Script did not complete. Exit all analysis."
  return 1
fi


# TODO
# Convert .h5ad into rds
# Substitute the file paths
. h5ad_to_rds.sh analysis/output/0_checkpoint/adata_final_Neuron_cca_features.h5ad analysis/output/0_checkpoint/adata_final_Neuron_cca_features.rds



# TODO
# TODO: Harmonize file name
FILE="S2.3_singlet_clustree_analysis.r"
Rscript "$FILE"
if [[ $? -eq 0 ]]
then
  echo "Success '$FILE'"
  echo "Script ran successfuly."
  rm -v ./output/0_checkpoint/adata_final_Neuron_cca_features.rds
else
  echo "Error '$FILE'"
  echo "Script did not complete. Exit all analysis."
  return 1
fi


# TODO
# TODO: Harmonize file name. Change file name to remove DR
# TODO: make it general for all populations
FILE="Cluster_fusion_DR.py"
python "$FILE"
if [[ $? -eq 0 ]]
then
  echo "Success '$FILE'"
  echo "Script ran successfuly."
else
  echo "Error '$FILE'"
  echo "Script did not complete. Exit all analysis."
  return 1
fi



# TODO
# TODO: file to visualize cluster fusion in UMAP

