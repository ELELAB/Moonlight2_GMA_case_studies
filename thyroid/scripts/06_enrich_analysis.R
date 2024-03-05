## This script performs enrichment analysis of three sets of genes:
## oncogenic mediators, functional gene-CpG pairs, and driver genes

## Load libraries ------------------------------------
library(Moonlight2R)
library(enrichR)
library(tidyverse)
library(patchwork)

## Source functions script  --------------------------
source("99_functions.R")

## Load data ------------------------------------

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
gene_lists <- list(oncogenic_mediators, driver_genes, functional_pairs_genes)
names(gene_lists) <- c("oncogenic_mediators", "driver_genes", "functional_pairs_genes")


## Analyze data ----------------------------------

## Perform enrichment analysis

# Create vector of databases to use for enrichment analysis 
dbs <- c("GO_Molecular_Function_2021",
         "GO_Biological_Process_2021",
         "MSigDB_Hallmark_2020")

# Perform enrichment analysis on each of the three gene lists
enrich_analysis <- map(gene_lists, function(x) { 
  enrichr(genes = x, databases = dbs)
})


## Visualize data ----------------------------------

# Visualize enrichment analysis results of GO molecular function terms of oncogenic mediators
mf_om <- goplot(data = enrich_analysis$oncogenic_mediators$GO_Molecular_Function_2021, 
                title = "GO Molecular Function 2021 - oncogenic mediators", 
                top = 10)

# Visualize enrichment analysis results of GO biological process terms of oncogenic mediators
bp_om <- goplot(data = enrich_analysis$oncogenic_mediators$GO_Biological_Process_2021, 
                title = "GO Biological Process 2021 - oncogenic mediators", 
                top = 10)

# Visualize enrichment analysis results of MSigDB Hallmark terms of oncogenic mediators
hm_om <- goplot(data = enrich_analysis$oncogenic_mediators$MSigDB_Hallmark_2020, 
                title = "MSigDB Hallmark 2020 - oncogenic mediators", 
                top = 10)

# Visualize enrichment analysis results of GO molecular function terms of driver genes
mf_dg <- goplot(data = enrich_analysis$driver_genes$GO_Molecular_Function_2021, 
                title = "GO Molecular Function 2021 - driver genes", 
                top = 10)

# Visualize enrichment analysis results of GO biological process terms of driver genes
bp_dg <- goplot(data = enrich_analysis$driver_genes$GO_Biological_Process_2021, 
                title = "GO Biological Process 2021 - driver genes", 
                top = 10)

# Visualize enrichment analysis results of MSigDB Hallmark terms of driver genes
hm_dg <- goplot(data = enrich_analysis$driver_genes$MSigDB_Hallmark_2020, 
                title = "MSigDB Hallmark 2020 - driver genes", 
                top = 10)

# Visualize enrichment analysis results of GO molecular function terms of functional pairs genes
mf_fp <- goplot(data = enrich_analysis$functional_pairs_genes$GO_Molecular_Function_2021, 
                title = "GO Molecular Function 2021 - functional pairs genes", 
                top = 10)

# Visualize enrichment analysis results of GO biological process terms of functional pairs genes
bp_fp <- goplot(data = enrich_analysis$functional_pairs_genes$GO_Biological_Process_2021, 
                title = "GO Biological Process 2021 - functional pairs genes", 
                top = 10)

# Visualize enrichment analysis results of MSigDB Hallmark terms of functional pairs genes
hm_fp <- goplot(data = enrich_analysis$functional_pairs_genes$MSigDB_Hallmark_2020, 
                title = "MSigDB Hallmark 2020 - functional pairs genes", 
                top = 10)

# Visualize enrichment analysis of each database of all gene lists side by side
mf_all <- mf_om + mf_fp 
bp_all <- bp_om + bp_dg + bp_fp
hm_all <- hm_om + hm_dg + hm_fp

# Visualize enrichment analysis of each database of Moonlight gene lists side by side
mf_moonlight <- mf_om
bp_moonlight <- bp_om + bp_dg
hm_moonlight <- hm_om + hm_dg


## Save data ------------------------------------

# Save results of enrichment analyses
save(enrich_analysis, file = "thyroid/results/06_enrich_analyses.rda")

# Save side by side enrichment plots of GO molecular function terms 
ggsave(filename = "06_enrich_GO_molecular_function.pdf",
       plot = mf_all,
       path = "thyroid/results/",
       width = 20,
       height = 10)

# Save side by side enrichment plots of GO biological process terms 
ggsave(filename = "06_enrich_GO_biological_process.pdf",
       plot = bp_all,
       path = "thyroid/results/",
       width = 20,
       height = 10)

# Save side by side enrichment plots of MSigDB hallmark terms 
ggsave(filename = "06_enrich_MSigDB_hallmark.pdf",
       plot = hm_all,
       path = "thyroid/results/",
       width = 20,
       height = 10)

# Save Moonlight side by side enrichment plots of GO molecular function terms 
ggsave(filename = "06_enrich_GO_molecular_function_moonlight.pdf",
       plot = mf_moonlight,
       path = "thyroid/results/",
       width = 12,
       height = 8)

# Save Moonlight side by side enrichment plots of GO biological process terms 
ggsave(filename = "06_enrich_GO_biological_process_moonlight.pdf",
       plot = bp_moonlight,
       path = "thyroid/results/",
       width = 12,
       height = 8)

# Save Moonlight side by side enrichment plots of MSigDB hallmark terms 
ggsave(filename = "06_enrich_MSigDB_hallmark_moonlight.pdf",
       plot = hm_moonlight,
       path = "thyroid/results/",
       width = 12,
       height = 8)




