## This script performs enrichment analysis of three sets of genes:
## oncogenic mediators, functional gene-CpG pairs, and driver genes

## Load libraries ------------------------------------
library(tidyverse)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(survMisc)

# Source functions
source("99_functions.R")


## Load data ------------------------------------

# Load PRA data (Moonlight's primary layer)
data_PRA <- get(load("lung/data/LUAD_dataPRA.rda"))

# Load functional gene-CpG pairs (EpiMix)
functional_pairs <- read_csv("lung/results/FunctionalPairs_Regular.csv")

# Load driver genes (Moonlight + EpiMix i.e. Moonlight's secondary layer)
data_EDA <- get(load("lung/results/LUAD_dataEDA.rda"))

# Load gene expression data
data_exp <- get(load("lung/data/LUAD_dataFilt_wrangled.rda"))

# Load clinical data
clinical_data <- get(load("lung/data/LUAD_clinical_data.rda"))

#clinical_data <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical") %>% 
#  as_tibble()


## Wrangle data ----------------------------------

# Get oncogenic mediators - putative TSGs and putative OCGs
oncogenic_mediators_TSGs <- data_PRA$TSG %>% names
oncogenic_mediators_OCGs <- data_PRA$OCG %>% names

# Get driver genes - TSGs and OCGs
driver_genes_TSGs <- data_EDA$TSG
driver_genes_OCGs <- data_EDA$OCG

# Get functional genes from EpiMix which are those functional genes that have 
# the same methylation state in all of its associated CpG probes
# and moreover all dual states are removed 
functional_pairs_genes <- functional_pairs %>%
  group_by(Gene) %>%
  dplyr::filter(all(State == first(State))) %>%
  ungroup() %>% 
  dplyr::filter(State != "Dual") 

# The putative TSGs among functional genes are those that have hypermethylated CpGs
# The putative OCGs among functional genes are those that have hypomethylated CpGs
functional_genes_TSGs <- functional_pairs_genes %>%
  group_by(Gene) %>% 
  dplyr::count(State) %>% 
  ungroup() %>% 
  dplyr::select(-c(n)) %>% 
  dplyr::filter(State == "Hyper") %>% 
  pull(Gene)
functional_genes_OCGs <- functional_pairs_genes %>%
  group_by(Gene) %>% 
  dplyr::count(State) %>% 
  ungroup() %>% 
  dplyr::select(-c(n)) %>% 
  dplyr::filter(State == "Hypo") %>% 
  pull(Gene)

# Merge all gene sets into one list
gene_sets <- list(oncogenic_mediators_TSGs, oncogenic_mediators_OCGs,
                  driver_genes_TSGs, driver_genes_OCGs,
                  functional_genes_TSGs, functional_genes_OCGs)
names(gene_sets) <- c("oncogenic_mediators_TSGs", "oncogenic_mediators_OCGs",
                      "driver_genes_TSGs", "driver_genes_OCGs",
                      "functional_genes_TSGs", "functional_genes_OCGs")

# Get vector of tumor samples from expression matrix
tumor_samples <- TCGAquery_SampleTypes(barcode = colnames(data_exp), 
                                       typesample = "TP") 

# Subset gene expression matrix to contain only tumor samples
data_exp_tumor <- data_exp[, tumor_samples]

# Convert column names of gene expression matrix to go from tumor sample
# level to patient level
colnames(data_exp_tumor) <- colnames(data_exp_tumor) %>% 
  str_extract(string = ., pattern = "^[^-]+-[^-]+-[^-]+")

# Subset the clinical data to contain only tumor patients of interest
# and select only variables related to survival
clinical_data_tumor <- clinical_data %>% 
  dplyr::filter(submitter_id %in% colnames(data_exp_tumor)) %>% 
  dplyr::select(submitter_id, days_to_last_follow_up,
                days_to_death, vital_status, gender, 
                age_at_diagnosis, ajcc_pathologic_stage) %>% 
  dplyr::mutate(time = case_when(vital_status == "Dead" ~ days_to_death,
                                 TRUE ~ days_to_last_follow_up),
                time_years = time/365,
                age_years = age_at_diagnosis/365,
                vital_status_binary = case_when(vital_status == "Dead" ~ 1,
                                                TRUE ~ 0),
                gender_binary = case_when(gender == "female" ~ 1,
                                          gender == "male" ~ 0),
                stage = 
                  case_when(str_detect(ajcc_pathologic_stage, "[ABC]$") ~ str_remove(ajcc_pathologic_stage, ".$"),
                            TRUE ~ ajcc_pathologic_stage),
                stage = case_when(stage == "Stage I" ~ 1,
                                  stage == "Stage II" ~ 2,
                                  stage == "Stage III" ~ 3,
                                  stage == "Stage IV" ~ 4))

# Apply a z-scale normalization on expression data 
data_exp_tumor_z <- t(data_exp_tumor) %>% 
  scale(., center = TRUE, scale = TRUE) %>% 
  t()



## Analyze data ---------------------------------- 

## Survival analysis 

## Test for Cox proportionality assumption before doing Cox regression analysis

# For each gene set i
cox_ph_test <- map(gene_sets, function(i) {
  
  # For each gene x in a gene set i,
  # test the Cox proportionality assumption
  map(i, function(x) {
    cox_reg_ph(exp_data = data_exp_tumor_z, gene = x, clinical_data = clinical_data_tumor)
  }) %>% 
    bind_rows()
  
})

## Proceed with only genes where the proportional hazards assumption is met

# For each tibble containing proportional hazards assumption test of a gene set
sig_genes <- map(cox_ph_test, function(x) {
  
  # Keep only genes that meet the assumption
  x %>% 
    group_by(gene) %>% 
    dplyr::filter(all(p > 0.05)) %>% 
    ungroup %>% 
    pull(gene) %>% 
    unique()
  
})


## Fit a univariate Cox regression model for each gene

# For each gene set i
cox_reg_uni_results <- map(sig_genes, function(i) {
  
  # For each gene x in a gene set i,
  # fit a univariate Cox regression model 
  cox_reg_gene_set <- map(i, function(x) {
    cox_reg_univariate(exp_data = data_exp_tumor_z, gene = x, clinical_data = clinical_data_tumor)
  }) %>% 
    bind_rows()
  
  # Correct p-values from univariate Cox regression for multiple testing using fdr method
  cox_reg_gene_set <- cox_reg_gene_set %>% 
    dplyr::mutate(fdr = p.adjust(`Pr(>|z|)`, method = "fdr"))
  
})


## Find genes with a significant effect on survival from univariate Cox
## regression and do multivariate Cox regression on those genes
sig_genes_multivariate <- map(cox_reg_uni_results, function(x) {
  
  # Get significant genes
  genes_for_multi <- x %>% 
    dplyr::filter(fdr < 0.05) %>% 
    pull(gene)
  
  # Perform multivariate Cox regression
  cox_multi <- map(genes_for_multi, function(i) {
    cox_reg_multivariate(exp_data = data_exp_tumor_z, gene = i, clinical_data = clinical_data_tumor) 
  }) %>% 
    bind_rows()
  
})

## Perform Kaplan-Meier survival analysis on genes with a significant
## effect on survival at multivariate level
km_sig_genes_multi <- map(sig_genes_multivariate, function(x) {
  
  # Extract genes in gene set x whose expression has a significant effect
  # on survival (at multivariate level)
  if (length(x) != 0) {
    sig_genes_exp <- x %>% 
      dplyr::filter(variable == "exp_value" & `Pr(>|z|)` < 0.05) %>% 
      pull(gene) %>% 
      unique  
    
    # Perform Kaplan-Meier analysis on above genes 
    km_results <- map(sig_genes_exp, function(i) {
      km_survival(exp_data = data_exp_tumor_z, gene = i, clinical_data = clinical_data_tumor)
    })
    names(km_results) <- sig_genes_exp
    return(km_results)
  }
})


## Overview of number of genes in each step of Cox regression analysis

# Get the number of genes in each gene set that had a significant effect on
# survival at the univariate level
num_genes_sig_uni <- map_int(sig_genes_multivariate, function(x) { if ("gene" %in% colnames(x)) {
  num_genes <- x %>% 
    dplyr::select(gene) %>% 
    pull %>% 
    unique %>% 
    length
  return(num_genes)
} else { 
  return(0) 
}}) %>% 
  as_tibble() %>% 
  pull

# Get the number of genes in each gene set that had a significant effect on
# survival at the multivariate level
num_genes_sig_multi <- map_int(sig_genes_multivariate, function(x) { if ("Pr(>|z|)" %in% colnames(x)) {
  num_genes <- x %>% 
    dplyr::filter(variable == "exp_value" & `Pr(>|z|)` < 0.05) %>% 
    pull %>% 
    unique %>% 
    length
  return(num_genes)
} else {
  return(0)
}}) %>% 
  as_tibble() %>% 
  pull

# Create overview of number of genes in each step of Cox regression analysis
num_genes <- tibble("gene_set" = names(gene_sets),
                    "total_num_genes" = map_dbl(gene_sets, length),
                    "total_num_genes_meet_assumption" = map_dbl(sig_genes, length),
                    "total_num_genes_sig_uni" = num_genes_sig_uni,
                    "total_num_genes_sig_multi" = num_genes_sig_multi) %>% 
  dplyr::mutate("percentage_sig_uni" = (total_num_genes_sig_uni / total_num_genes) * 100,
                "percentage_sig_multi" = (total_num_genes_sig_multi / total_num_genes) * 100)


## Save data ------------------------------------

# Save results of proportional hazards assumption test
save(cox_ph_test, file = "lung/results/09_LUAD_cox_ph.rda")

# Save results of Cox univariate regression analysis
save(cox_reg_uni_results, file = "lung/results/09_LUAD_cox_univariate.rda")

# Save results of Cox multivariate regression analysis
save(sig_genes_multivariate, file = "lung/results/09_LUAD_cox_multivariate.rda")

# Save overview of number of genes in each step of Cox regression
write_csv(num_genes, file = "lung/results/09_LUAD_cox_num_genes.csv")

# Save all Kaplan-Meier survival plots of driver genes
km_sig_genes_multi_drivers <- km_sig_genes_multi[c("driver_genes_OCGs",
                                                   "driver_genes_TSGs")]
map(km_sig_genes_multi_drivers, function(x) {
  
  if (!is.null(x) & length(x) != 0) {
    map(seq.int(length(x)), function(i) {
      
      pdf(str_c("lung/results/09_survival_analysis_plots/09_LUAD_km_survival_", names(x)[i], ".pdf"),
          width = 10, height = 10)
      print(x[[i]], newpage = FALSE)
      dev.off()
      
    })
  }
})




