## This script filters for patient barcodes containing methylation, but no mutations

## Load libraries
library(tidyverse)
library(stringr)

## Load data------------
load("breast_basal/data/BRCA_basal_mutation_subset.rda")
load("breast_basal/results/DEG_Mutations_Annotations.rda")
load("breast_basal/results/Oncogenic_mediators_mutation_summary.rda")

## Fetch names of driver genes predicted by DMA 
DMA_Driver_genes <- Oncogenic_mediators_mutation_summary %>% 
  filter(CScape_Driver >= 1 ) %>% 
  select(Hugo_Symbol)

## Fetch patient barcodes containing these driver genes from the GMA dataset
DMA_within_GMA <- DEG_Mutations_Annotations %>% 
  filter(Hugo_Symbol %in% DMA_Driver_genes$Hugo_Symbol) %>% 
  mutate(Tumor_Sample_Barcode_short = substr(x = Tumor_Sample_Barcode, start = 1, stop = 15))

## Check if any patient sample contains driver genes predicted by GMA 
## but is missing from the DMA_within_GMA barcodes 
Missing_sample_in_DMA <-all_mutations_subset %>%
  filter(Tumor_Sample_Barcode %in% DMA_within_GMA$Tumor_Sample_Barcode_short)
