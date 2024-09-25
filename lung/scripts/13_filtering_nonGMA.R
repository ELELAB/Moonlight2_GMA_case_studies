## This script filters for patient barcodes containing methylation, but no mutations in the MAF file

## Load libraries
library(tidyverse)
library(stringr)

## Load data------------
load("lung/data/LUAD_mutations_subset.rda")
load("lung/data/LUAD_sample_annotations.rda")

##Load patient barcodes having mutations
patients_mutated <- mutations_subset %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  select(Tumor_Sample_Barcode) %>% 
  unique() %>% 
  mutate(primary = substr(x = Tumor_Sample_Barcode, start = 1, stop = 15)) %>% 
  ungroup %>% 
  select(primary)

#Load patient barcodes having methylations
patients_methylated <- sample_annotation %>% 
  select(primary)

#Find patients having methylations but no mutations
methylated_but_not_mutated <- setdiff(patients_methylated, patients_mutated)
