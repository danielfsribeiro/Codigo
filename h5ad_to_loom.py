import sys
import scanpy as sc


if len(sys.argv) < 3:
    print("Usage:")
    print("python3 h5ad_to_loom.py path/to/h5ad path/to/loom")
    exit(0)

source_file = sys.argv[1]
dest_file = sys.argv[2]
print(source_file)
print(dest_file)

print("Load h5ad...")
adata = sc.read_h5ad(source_file)
if source_file.find('seurat') != -1:
    adata.var_names.name = None
    adata.obs_names.name = None
print("Write loom...")
adata.write_loom(dest_file, write_obsm_varm=False)  # no obsm will be saved
