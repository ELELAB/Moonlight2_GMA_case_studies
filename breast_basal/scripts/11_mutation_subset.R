## This script subsets the patient barcodes with mutations from the entire methylation dataset

## Load libraries-----------------
library(tidyverse)
library(dplyr)
library(stringr)


## Load data-----------
Maf_Basal <- read.csv("breast_basal/data/mutations.csv")
sample_annotation <- get(load("breast_basal/data/BRCA_Basal_sample_annotations.rda"))

## DMA----------
Maf_Basal_filtered <- Maf_Basal |> 
  mutate(VAF = t_alt_count/t_depth) %>% 
  filter(VAF > 0.05, 
         t_alt_count >= 3,
         t_depth > 30,
         n_depth >10)

## Extract the cancer samples
cancer_barcodes <- sample_annotation  |> 
  filter(sample.type == "Cancer")  |> 
  select(primary)

## Match methylation patient barcodes with mutations for the cancer samples
mutations_subset <-  Maf_Basal_filtered %>% 
  filter((substr(x = Tumor_Sample_Barcode, start = 1, stop = 15)) %in% cancer_barcodes$primary)

## Save the mutation subset
save(mutations_subset,file='breast_basal/data/BRCA_basal_mutation_subset.rda')
