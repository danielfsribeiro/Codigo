#!/bin/bash
if [[ ${#@} -ne 2 ]]; then
    echo "Usage:"
    echo ". h5ad_to_rdata.sh path/to/h5ad path/to/rdata"
    return 0
fi

h5ad_file=$1
shift
rdata_file=$1

# Check files and folders exist
if [[ -e $h5ad_file ]] && [[ -d ${rdata_file%/*} ]]; then
    conda activate py39
    python h5ad_to_loom.py $h5ad_file "${h5ad_file%.h5ad}.loom"
    Rscript loom_to_rds.r "${h5ad_file%.h5ad}.loom" $rdata_file
    rm -v "${h5ad_file%.h5ad}.loom"
else
    echo "File '$h5ad_file' or directory '${rdata_file%/*}' are not available"
fi