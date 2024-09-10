## This script runs the DMA function on the mutation subset of the methylation samples

## Load R environment-------------
source('renv/activate.R')

## Load libraries----------------
library(Moonlight2R) 
library(tidyverse)

## Load the mutation subset along with the DEA and PRA objects from the GMA analysis---------------------
DEA_Basal <- get(load("breast_basal/data/BRCA_basal_sig_DEA.rda"))
mutations_subset <- get(load("breast_basal/data/BRCA_basal_mutation_subset.rda"))
PRA_Basal <- get(load("breast_basal/data/BRCA_basal_dataPRA.rda"))

## Run the DMA function
DMA_Basal <- DMA(dataMAF = mutations_subset,
                dataDEGs = DEA_Basal,
                dataPRA = PRA_Basal,
                results_folder = "breast_basal/results/",
                coding_file = Sys.getenv('CSCAPE_CODING'),
                noncoding_file = Sys.getenv('CSCAPE_NONCODING'))

## Save the DMA data object
save(DMA_Basal, file = "breast_basal/data/BRCA_basal_dataDMA.rda")

