## This script subsets the patient barcodes with mutations from the entire methylation dataset

## Load libraries-----------------
library(tidyverse)
library(stringr)

## Load data-----------
Maf_Basal <- read.csv("lung/data/mutations.csv")
load('lung/data/LUAD_sample_annotations.rda')

## DMA----------
Maf_Basal_filtered <- Maf_Basal |> 
  mutate(VAF = t_alt_count/t_depth) %>% 
  filter(VAF > 0.05, 
         t_alt_count >= 3,
         t_depth > 30,
         n_depth >10)

## Separate the samples into cancer and normal-------
# cancer_barcodes <- sample_annotation  |> 
#  filter(sample.type == "Cancer")  |> 
#  select(primary)
# normal_barcodes <- sample_annotation  |> 
#  filter(sample.type == "Normal")  |> 
#  select(primary)

## Match methylation patient barcodes with mutations for the cancer samples--------
all_mutations_subset <-  Maf_Basal_filtered %>% 
  filter((substr(x = Tumor_Sample_Barcode, start = 1, stop = 15)) %in% sample_annotation$primary)

## Match methylation patient barcodes with mutations for the normal samples----------
# normal_mutations_subset <- Maf_Basal_filtered %>% 
#  filter((substr(x = Matched_Norm_Sample_Barcode, start = 1, stop = 15)) %in% normal_barcodes$primary)

## Combine the barcode-matched cancer and normal samples---------
#all_mutations_subset <- cancer_mutations_subset %>% 
#  union(y = normal_mutations_subset)

## Save the combined mutation subset-------
save(all_mutations_subset, file = "lung/data/LUAD_mutations_subset.rda")
