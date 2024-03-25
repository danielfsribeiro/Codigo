import scanpy as sc

Astrocytes = "/mnt/c/Users/vasco/OneDrive/Área de Trabalho/Bioinformatica_IMM/2023-07-16_1809/Astrocytes_imune_Cells/adata_final_Astrocyte_cca_features.h5ad"

Immune = "/mnt/c/Users/vasco/OneDrive/Área de Trabalho/Bioinformatica_IMM/2023-07-16_1809/Astrocytes_imune_Cells/adata_final_Immune_cca_features.h5ad"

Astrocytes_adata = sc.read(Astrocytes)
Immune_adata = sc.read(Immune)

print("Contents of Astrocytes adata:")
print(Astrocytes_adata)

print("\nContents of Immune adata:")
print(Immune_adata)

