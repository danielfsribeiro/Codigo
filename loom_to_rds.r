library("Seurat")
library("SeuratDisk")
library("SeuratData")


args = commandArgs(trailingOnly = TRUE)
print(args)
if (length(args) < 2)
{
    print("Usage:")
    print("Rscript loom_to_rdata.r path/to/loom path/to/Rdata")
    quit(save = 'no')
}

source_file = args[1]
dest_file = args[2]

print("Load loom...")
adata = Connect(source_file)
print("Convert loom to Seurat...")
adata.seurat = as.Seurat(adata, features = "var_names", cells = "obs_names")
saveRDS(adata.seurat, file = dest_file)

adata$close_all()

print("*** DONE ***")