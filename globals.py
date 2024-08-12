# Global definitions imported by multiple analysis files

H5AD_DIR = "/mnt/c/Users/vasco/Env/Astrocytes_imune_Cells"
CLUSTER_FUSION_DIR = "/mnt/c/Users/vasco/Env/cluster_fusion"
DATASETS = ["Immune", "Astrocyte"]
Clustresults = "/mnt/c/Users/vasco/Env/Astrocytes_imune_Cells/Clustresults"
ENV = "/mnt/c/Users/vasco/Env"


#H5AD_DIR = "Astrocytes_imune_Cells"
#CLUSTER_FUSION_DIR = "cluster_fusion"
#DATASETS = ["Immune", "Astrocyte"]
#Clustresults = "Astrocytes_imune_Cells/Clustresults"
#ENV = "./"

# Lineage resolution tests for clustree
lineage_resolution_tests = {
    'Astrocyte': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1, 2],
    'Immune': [0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.1, 1.2, 2]}
