## This script runs GMA on gene expression and methylation data

## Load libraries ------------------------------------
library(EpiMix)
library(Moonlight2R)
library(tibble)
library(tidyr)
library(readr)
library(fuzzyjoin)
library(stringr)
library(magrittr)
library(BiocGenerics)
library(sesameData)


## Load data ------------------------------------

# Load filtered gene expression data subsetted to contain same patients as 
# methylation data
data_exp <- get(load("breast_basal/data/BRCA_basal_dataFilt_wrangled.rda"))

# Load table of DEGs
data_DEGs <- get(load("breast_basal/data/BRCA_basal_sig_DEA.rda"))

# Load methylation data
data_met <- get(load("breast_basal/data/BRCA_basal_dataMet_wrangled.rda")) 

# Load sample annotation
sample_annotations <- get(load("breast_basal/data/BRCA_Basal_sample_annotations.rda"))

# Load oncogenic mediators predicted from Moonlight's primary layer
data_PRA <- get(load("breast_basal/data/BRCA_basal_dataPRA.rda"))

# Load data from Moonlight package
data(MetEvidenceDriver)
data(NCG)
data(EncodePromoters)


## Analyze data ------------------------------------

# Run GMA on data
data_EDA <- GMA(dataMET = data_met, 
                dataEXP = data_exp, 
                dataPRA = data_PRA, 
                dataDEGs = data_DEGs, 
                sample_info = sample_annotations, 
                met_platform = "HM450", 
                prevalence_filter = NULL, 
                output_dir = "breast_basal/results/", 
                cores = 1, 
                roadmap.epigenome.ids = "E119", 
                roadmap.epigenome.groups = NULL)


# Save data ------------------------------------

# Save driver genes predicted using EDA
save(data_EDA, 
     file = "breast_basal/results/BRCA_basal_dataEDA.rda")
