#!/bin/bash
# This script organizes the multiple analysis files

conda activate py39

H5AD_DIR="/mnt/c/Users/vasco/Env/Astrocytes_imune_Cells"
#H5AD_DIR="Astrocytes_imune_Cells"


## Creates cluslters at multiple resolutions
#FILE="Leiden script.py"
#python "$FILE"
#if [[ $? -eq 0 ]]
#then
#  echo "Success '$FILE'"
#  echo "Script ran successfuly."
#  echo " "
#else
#  echo "Error '$FILE'"
#  echo "Script did not complete. Exit all analysis."
#  echo " "
#  return 1
#fi

# Plot UMAPs
#FILE="Astrocytes_GoodResolution.py"
#python "$FILE"
#if [[ $? -eq 0 ]]
#then
#  echo "Success '$FILE'"
#  echo "Script ran successfuly."
#  echo " "
#else
#  echo "Error '$FILE'"
#  echo "Script did not complete. Exit all analysis."
#  echo " "
#  return 1
#fi

#FILE="Immune_GoodResolution.py"
#python "$FILE"
#if [[ $? -eq 0 ]]
#then
#  echo "Success '$FILE'"
#  echo "Script ran successfuly."
#  echo " "
#else
#  echo "Error '$FILE'"
#  echo "Script did not complete. Exit all analysis."
#  echo " "
#  return 1
#fi

# Convert .h5ad into rds
#. h5ad_to_rds.sh ${H5AD_DIR}/adata_final_Astrocyte_cca_features.h5ad ${H5AD_DIR}/adata_final_Astrocyte_cca_features.rds
#. h5ad_to_rds.sh ${H5AD_DIR}/adata_final_Immune_cca_features.h5ad ${H5AD_DIR}/adata_final_Immune_cca_features.rds
#FILE="Clustree_analysis.r"
#Rscript "$FILE"
#if [[ $? -eq 0 ]]
#then
#  echo "Success '$FILE'"
#  echo "Script ran successfuly."
#  echo " "
#  rm -v ${H5AD_DIR}/adata_final_Astrocyte_cca_features.rds
#  rm -v ${H5AD_DIR}/adata_final_Immune_cca_features.rds
#else
#  echo "Error '$FILE'"
#  echo "Script did not complete. Exit all analysis."
#  echo " "
#  return 1
#fi


# TODO
# TODO: Harmonize file name. Change file name to remove DR
# TODO: make it general for all populations
FILE="Cluster_fusion.py"
python "$FILE"
if [[ $? -eq 0 ]]
then
  echo "Success '$FILE'"
  echo "Script ran successfuly."
  echo " "
else
  echo "Error '$FILE'"
  echo "Script did not complete. Exit all analysis."
  echo " "
  return 1
fi



# TODO
# TODO: file to visualize cluster fusion in UMAP

