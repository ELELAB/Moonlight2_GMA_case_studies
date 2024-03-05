## This script compares EpiMix and Moonlight results with the Cancer Gene Census

## Load libraries ------------------------------------
library(tidyverse)
library(viridis)


## Load data ------------------------------------

# Load gene expression matrix
data_exp <- get(load("thyroid/data/THCA_dataFilt_wrangled.rda"))

# Load cancer gene census data
cgc <- read_csv("thyroid/data/cancer_gene_census.csv")

# Load PRA data (Moonlight's primary layer)
data_PRA <- get(load("thyroid/data/THCA_dataPRA.rda"))

# Load functional gene-CpG pairs (EpiMix)
functional_pairs <- read_csv("thyroid/results/FunctionalPairs_Regular.csv")

# Load driver genes (Moonlight + EpiMix i.e. Moonlight's secondary layer)
data_EDA <- get(load("thyroid/results/THCA_dataEDA.rda"))


## Wrangle data ----------------------------------

# Get oncogenic mediators
oncogenic_mediators <- Reduce(c, data_PRA) %>% names

# Get driver genes
driver_genes <- Reduce(c, data_EDA)

# Get genes from EpiMix which are those functional genes that have 
# the same methylation state in all of its associated CpG probes
# and moreover all dual states are removed 
functional_pairs_genes <- functional_pairs %>%
  group_by(Gene) %>%
  dplyr::filter(all(State == first(State))) %>%
  ungroup() %>% 
  dplyr::filter(State != "Dual") %>% 
  dplyr::select(Gene) %>% 
  pull %>% 
  unique

# Wrangle the three gene lists into one list
gene_lists <- list(oncogenic_mediators, functional_pairs_genes, driver_genes)
names(gene_lists) <- c("Moonlight's primary layer", "EpiMix", "EDA")


## Analyze data ------------------------------------

# Compute various quantitative statistics to compare gene sets with CGC
overlaps <- map(seq.int(length(gene_lists)), function(x) {
  
  # Calculate true positives which is the overlap between CGC and gene list x (absolute number)
  TP <- intersect(gene_lists[[x]], cgc$`Gene Symbol`) %>% length
  
  # Calculate false positives which is genes in gene set not identified in CGC (absolute number)
  FP <- setdiff(gene_lists[[x]], cgc$`Gene Symbol`) %>% length
  
  # Calculate false negatives which is genes in CGC but not identified in gene set (absolute number)
  FN <- setdiff(cgc$`Gene Symbol`, gene_lists[[x]]) %>% length
  
  # Calculate genes that are not in gene set and not in CGC
  # To be used in Fisher's exact test
  # N is estimate for the total number of genes which would be
  # the union of genes in CGC and genes in exp matrix
  N <- union(cgc$`Gene Symbol`, rownames(data_exp)) %>% length
  non_driver_cgc <- N - length(union(gene_lists[[x]], cgc$`Gene Symbol`))
  
  # Calculate precision which is TP/(TP+FP) * 100 (percentage)
  precision <- TP / (TP + FP) * 100
  
  # Calculate sensitivity which is TP/(TP+FN) * 100 (percentage)
  sensitivity <- TP / (TP + FN) * 100
  
  # Perform Fisher's exact test and extract results from this test
  mat_for_fisher_test <- matrix(c(TP, FP, FN, non_driver_cgc), nrow = 2)
  ft <- fisher.test(mat_for_fisher_test)
  ft_pvalue <- ft$p.value
  ft_odds_ratio <- ft$estimate
  
  # Combine results in table
  overlap_tbl <- tibble(gene_set = names(gene_lists)[[x]], 
                        gene_set_num = length(gene_lists[[x]]),
                        TP = TP,
                        FP = FP,
                        FN = FN,
                        precision = precision,
                        sensitivity = sensitivity,
                        fisher_pvalue = ft_pvalue,
                        fisher_odds_ratio = ft_odds_ratio) 
}) %>% 
  bind_rows() 


## Save data ------------------------------------

# Save table showing gene overlap numbers
write_csv(overlaps, file = "thyroid/results/07_gene_overlap_cgc.csv")






