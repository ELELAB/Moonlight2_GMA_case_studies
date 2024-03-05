
Cancer Systems Biology, Section of Bioinformatics, Department of Health Technology, Technical University of Denmark, 2800, Lyngby, Denmark

# Moonlight2_GMA_case_studies

This repository contains scripts for case studies related to the discovery of cancer driver genes using the Moonlight
framework. The case studies are conducted on basal-like breast cancer, lung adenocarcinoma, and thyroid cancer using 
data from The Cancer Genome Atlas (TCGA). The associated publication is:

bioRxiv reference

Please cite the above publication if you use the contents, scripts or results for your own research.

Below are instructions for reproducing the analyses. 

## Structure and content of GitHub and OSF repositories

This GitHub repository contains scripts associated with the publication
with a main folder for each cancer (sub)type. Within each cancer (sub)type folder,
a subfolder called `scripts` contains the associated scripts. The scripts are
numbered according to the order in which they are run.

The corresponding [OSF repository](https://osf.io/j4n8q/) contains data and results associated with
the scripts and is organized in the same way as the GitHub repository with a
main folder for each cancer (sub)type. Within each cancer (sub)type folder,
subfolders called `data` and `results` contain the associated data and results,
respectively. The results files in `results` are numbered according to the
script that generated them.

## Installing requirements and reproducing the analysis

All the analyses have been performed on a GNU/Linux server.

NB: When reproducing the analyses and results, the user cannot expect to obtain identical results to the ones
of the case studies and associated with the publication due to stochastic processes in the GRN step of the Moonlight protocol. 

### Computing environment

In order to reproduce the paper data, you will need to set up a `conda` environment
on which the expected version of R and the required packages will be installed;
this requires being able to run Anaconda by means of the `conda` executable.

If you don't have access to `conda` please see the [Miniconda installer page](https://docs.conda.io/en/latest/miniconda.html) for instructions on how to install Miniconda.

Once you have access to `conda`, follow the below instructions:

1. Clone our github repository into a local directory on your local machine:

```
git clone https://github.com/ELELAB/Moonlight2_GMA_case_studies.git
cd Moonlight2_GMA_case_studies
```

2. Create a virtual environment using conda and activate it. 
The environment directory should be placed in the Moonlight2_GMA_case_studies folder:

```
conda create --prefix ./methyl_case -c conda-forge r-base=4.3 r-pacman=0.5.1 r-curl=5.1.0 r-ragg=1.2.7 r-renv=1.0.5 r-osfr=0.2.9 r-devtools=2.4.5 gsl=2.7 gmp=6.2.1 glpk=5.0
conda activate ./methyl_case
```

3. Download data from the [COSMIC Cancer Gene Census](https://cancer.sanger.ac.uk/census).
This data can be downloaded from https://cancer.sanger.ac.uk/census by exporting it as
a `csv` file or from https://cancer.sanger.ac.uk/cosmic/download/cosmic/v99/cancergenecensus
by choosing the file from the `CRCh28` genome and afterwards converting it to a `csv` file. 
Once the data from the Cancer Gene Census has been downloaded, it must be a `csv` file named 
`cancer_gene_census.csv` and this file must be placed in the `data` folder of each cancer (sub)type:

```
breast_basal/data/cancer_gene_census.csv
lung/data/cancer_gene_census.csv
thyroid/data/cancer_gene_census.csv
```

4. Run the analyses:

```
bash ./run_all.sh
```

**WARNING**: our scripts use the [renv](https://rstudio.github.io/renv/articles/renv.html)
R package to handle automatic dependency installation. `Renv` writes packages in
its own cache folder, which is by default in the user's home directory. This might not be
desirable if free space in the home directory is limited. You can change the location of
the `Renv` root folder by setting a system environment variable - please see comments
in the `run_all.sh` script.

The `run_all.sh` script will perform the following steps to reproduce all results and data:

1. Download data from the corresponding [OSF repository](https://osf.io/j4n8q/) which
contains the required data to run the analyses and all results associated with the analyses.

2. Install in the environment all necessary packages to run the analyses.

4. Perform all analyses for basal-like breast cancer.

5. Perform all analyses for lung adenocarcinoma. 

6. Perform all analsyes for thyroid cancer. 

7. Compare results across cancer (sub)types. 
