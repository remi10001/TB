# This script makes a directory and downloads data from the ARCHS4 resource for the paper.
# run from the repository base directory using "bash download_ARCHS4_datatable.sh"

# Downloading this matrix table (4.0 GB) will take some time

mkdir data/ARCHS4
cd data/ARCHS4

wget https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5