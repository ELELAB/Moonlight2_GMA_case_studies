## This script filters for patient barcodes containing methylation, but no mutations

## Load libraries
library(tidyverse)
library(stringr)

## Load data------------
load("breast_basal/data/BRCA_basal_mutation_subset.rda")
load("breast_basal/data/BRCA_Basal_sample_annotations.rda")

## Fetch tumor sample barcodes from cancer mutation subset
tumor_patient_barcodes <- mutations_subset %>% 
  select(Tumor_Sample_Barcode) %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  unique() %>% 
  mutate(Tumor_Sample_Barcode_short = substr(x = Tumor_Sample_Barcode, start = 1, stop = 15))

## Check if any patient sample contains driver genes predicted by GMA 
## but is missing in MAF mutation file 
Missing_sample_in_MAF <- data.frame(setdiff(sample_annotation$primary,tumor_patient_barcodes$Tumor_Sample_Barcode_short))
