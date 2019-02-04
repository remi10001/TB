# TB

This repository provides code to reproduce the analyses of the unpublished manuscript titled "_Blood RNA Signatures Predict Recent Tuberculosis Exposure in Mice, Macaques and Humans_". The major analyses of the paper can be reproduced by running the jupyter notebook "Time Since Infection_Submission_Streamlinedcode.ipynb". All code written during the development of this project is in "Time Since Infection RNA_Paper_Analysis.ipynb" and may be reviewed. However, it is not structured for easy, streamlined running of all code.

This code was developed and used in a linux environment and has not been modified to aid in running with other platforms.

### This software is available under the MIT license. However, the script utils_immunoStates.R contains code that is copied from the MetaIntegrator R package, and thus this script is under an LGPL license.

## Installation (for running the provided jupyter notebook with R; alternatively, one can use the provided R markdown file with a local installation of R studio, or simply view the provided .html)

1. Install anaconda (a local installation manager and python distribution)[https://docs.anaconda.com/anaconda/install/], and create a new conda environment

```
$  conda create --name TB  python=3.6.3  
$  conda activate TB
```

2. Install or make sure the following python packages are installed:

- Pandas
- h5py
- jupyter
- jupyter notebook extensions (for easy browsing of the analysis notebook)

```
conda install pandas h5py jupyter
conda install -c conda-forge jupyter_contrib_nbextensions
```

3. Install r-essentials with mro base R version 3.4.3 (the R version used in this analysis)

```
conda install -c r mro-base=3.4.3 r-essentials 
```

4. Open R from command-line and run IRkernel so that jupyter notebook can see the R installation

```
R
> IRkernel::installspec()
```

5. Alternatively, use the provided env_list.txt and environment.yml files to reproduce the environment in which this analysis was performed. However, many packages in this environment are not necessary.

## Download required public datasets

1. Run the provided download scripts from the base repository directory to download required public data for the analyses.
```
bash download_macaque_human_data.sh
bash download_ARCHS4_datatable.sh
```

## Perform the analyses

1. Reproduce the major analyses of the paper by running the jupyter notebook "Time Since Infection_Submission_Streamlinedcode.ipynb" and following the steps it describes.