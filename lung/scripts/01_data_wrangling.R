## This script wrangles and prepares data for analysis

## Load libraries ------------------------------------
library(magrittr)
library(stringr)
library(TCGAbiolinks)
library(EpiMix)
library(tibble)
library(dplyr)

source("EpiMix_functions.R")

## Load data ------------------------------------

# Load filtered gene expression data
data_exp <- get(load("lung/data/LUAD_dataFilt_HUGO.rda"))

# Load methylation data
data_met <- get(load("lung/data/LUAD_DNAMet.rda"))


## Wrangle data ------------------------------------

# Subset columns (sample barcodes) to be on patient level
colnames(data_exp) <- colnames(data_exp) %>%
  str_extract("^([^\\-]+-[^\\-]+-[^\\-]+-[^\\-]+)") %>%
  str_remove(".$")
colnames(data_met) <- colnames(data_met) %>%
  str_replace_all("\\.", "-") %>% 
  str_extract("^([^\\-]+-[^\\-]+-[^\\-]+-[^\\-]+)") %>%
  str_remove(".$")

# Create annotation table that annotates each sample in methylation data
# as either cancer or normal
TP_samples <- TCGAquery_SampleTypes(colnames(data_met), 
                                    "TP")
NT_samples <- TCGAquery_SampleTypes(colnames(data_met), 
                                    "NT")
TP_samples_anno <- data.frame("primary" = TP_samples, 
                              "sample.type" = "Cancer") 
NT_samples_anno <- data.frame("primary" = NT_samples, 
                              "sample.type" = "Normal") 
sample_annotation_full <- rbind(TP_samples_anno, 
                                NT_samples_anno)

# Process methylation data
data_met_prep <- Preprocess_DNAMethylation(methylation.data = data_met, 
                                           met.platform = "HM450", 
                                           sample.info = sample_annotation_full, 
                                           group.1 = "Cancer", 
                                           group.2 = "Normal")

# Subset expression and methylation data matrices to contain the same patients
patients_intersection <- intersect(colnames(data_exp), 
                                   colnames(data_met_prep))
data_exp_subsetted <- data_exp[, patients_intersection]
data_met_subsetted <- data_met_prep[, patients_intersection] 

# Update annotation table that annotates each sample in methylation data
# as either cancer or normal
TP_samples <- TCGAquery_SampleTypes(colnames(data_met_subsetted), 
                                    "TP")
NT_samples <- TCGAquery_SampleTypes(colnames(data_met_subsetted), 
                                    "NT")
TP_samples_anno <- data.frame("primary" = TP_samples, 
                              "sample.type" = "Cancer") 
NT_samples_anno <- data.frame("primary" = NT_samples, 
                              "sample.type" = "Normal") 
sample_annotation <- rbind(TP_samples_anno, 
                           NT_samples_anno)


## Save data ------------------------------------

# Save wrangled filtered gene expression data
save(data_exp_subsetted,
     file = "lung/data/LUAD_dataFilt_wrangled.rda")

# Save wrangled methylation data
save(data_met_subsetted, 
     file = "lung/data/LUAD_dataMet_wrangled.rda")

# Save sample annotations
save(sample_annotation,
     file = "lung/data/LUAD_sample_annotations.rda")

