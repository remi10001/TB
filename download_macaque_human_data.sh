# This script makes directories and downloads data for the paper.
# run from the repository base directory using "bash download_macaque_human_data.sh"

mkdir data/GSE124688

mkdir data/GSE84152
cd data/GSE84152
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84152/suppl/GSE84152_non-normalized.txt.gz

gzip -d GSE84152_non-normalized.txt.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84152/matrix/GSE84152_series_matrix.txt.gz


cd ../..
mkdir data/GSE94438
cd data/GSE94438

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94438/matrix/GSE94438_series_matrix.txt.gz

wget  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94438/suppl/GSE94438_rawCounts_GeneNames_AllSamples.csv.gz
gzip -d GSE94438_rawCounts_GeneNames_AllSamples.csv.gz

cd ../..
mkdir data/GSE79362
cd data/GSE79362

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79362/matrix/GSE79362_series_matrix.txt.gz

cd ../..
mkdir data/GSE116014
cd data/GSE116014

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116014/matrix/GSE116014_series_matrix.txt.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116014/suppl/GSE116014_RAW.tar

tar xf GSE116014_RAW.tar
gzip -d *idat.gz