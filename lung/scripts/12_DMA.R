## This script runs the DMA function on the mutation subset of the methylation samples

## Load R environment-------------
source('renv/activate.R')

## Load libraries----------------
library(Moonlight2R) 
library(tidyverse)

## Load the mutation subset along with the DEA and PRA objects from the GMA analysis---------------------
DEA_lung <- get(load("lung/data/LUAD_sig_DEA.rda"))
mutations_subset <- get(load("lung/data/mutations_subset.rda"))
PRA_lung <- get(load("lung/data/LUAD_dataPRA.rda"))

## Run the DMA function---------
DMA_lung <- DMA(dataMAF = mutations_subset,
                 dataDEGs = DEA_lung,
                 dataPRA = PRA_lung,
                 results_folder = "lung/results/",
                 coding_file = Sys.getenv('CSCAPE_CODING'),
                 noncoding_file = Sys.getenv('CSCAPE_NONCODING'))

## Save the DMA data object----------
save(DMA_lung, file = "lung/data/LUAD_dataDMA.rda")




