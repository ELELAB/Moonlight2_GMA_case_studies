
Cancer Systems Biology, Section of Bioinformatics, Department of Health Technology, Technical University of Denmark, 2800, Lyngby, Denmark

# Moonlight2_GMA_case_studies

This repository contains scripts for case studies related to the discovery of cancer driver genes using the Moonlight
framework. The case studies are conducted on basal-like breast cancer, lung adenocarcinoma, and thyroid cancer using 
data from The Cancer Genome Atlas (TCGA). The associated publication is:

bioRxiv reference

Please cite the above publication if you use the contents, scripts or results for your own research.

Below are instructions on the first steps for reproducing the analyses. Afterwards, please see README files in each 
cancer folder for specifics on how to reproduce analyses for the given cancer type.

## Installing requirements and reproducing the analysis

All the analyses have been performed on a GNU/Linux server.

NB: When reproducing the analyses and results, the user cannot expect to obtain identical results to the ones
of the case studies and associated with the publication due to stochastic processes in the GRN step of the Moonlight protocol. 

### Computing environment

In order to reproduce the paper data, you will need to set up a `conda` environment
on which the expected version of R and the required packages will be installed;
this requires being able to run Anaconda by means of the `conda` executable.

If you don't have access to `conda` please see the [Miniconda installer page](https://docs.conda.io/en/latest/miniconda.html) for instructions on how to install Miniconda.

Once you have access to `conda`, you can

1. clone our github repository into a local directory on your local machine:

```
git clone https://github.com/ELELAB/Moonlight2_GMA_case_studies.git
cd Moonlight2_GMA_case_studies
```

2. create a virtual environment using conda and activate it. 
The environment directory should be placed in the Moonlight2_GMA_case_studies folder:

```
conda create --prefix ./methyl_case -c conda-forge r-base=4.3 r-pacman=0.5.1 r-curl=5.1.0 r-ragg=1.2.7 r-renv=1.0.3 r-osfr=0.2.9 r-devtools=2.4.5 gsl=2.7 gmp=6.2.1 glpk=5.0
conda activate ./methyl_case
```

3. change directory to the cancer type of interest, i.e. one of below possibilites:

```
cd breast_basal/scripts
cd lung/scripts
cd thyroid/scripts
```

4. follow instructions in the README of the specific cancer type to reproduce analyses
associated with that cancer type

