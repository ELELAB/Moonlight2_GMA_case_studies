## This script performs enrichment analysis on the drivers predicted by both GMA and DMA

## Load libraries ------------
library(tidyverse)
library(RColorBrewer)
library(enrichR)

## Source functions script -----------
source("99_functions.R")

## Load the GMA and DMA data  --------------------------
data_EDA <- get(load("breast_basal/results/BRCA_basal_dataEDA.rda"))
data_DMA <- get(load("breast_basal/data/BRCA_basal_DMA.rda"))

# Get driver genes
GMA_driver_genes <- Reduce(c, data_EDA)
DMA_driver_genes <- Reduce(c, data_DMA)

# Find overlapping genes
overlap <- intersect(GMA_driver_genes, DMA_driver_genes)

# Wrangle genes into a list 
overlap_list <- list(overlap)
names(overlap_list) <- c("Overlaps")

# Load MSigDB Hallmarks
dbs <- c("MSigDB_Hallmark_2020")

## Analyze data ------------

# Perform enrichment analysis
enrich_analysis <- map(overlap_list, function(x) { 
  enrichr(genes = x, databases = dbs)
})

# Plot enrichment analysis results
overlap_dr <- goplot(data = enrich_analysis$Overlaps$MSigDB_Hallmark_2020, 
                     title = "Overlap between DMA and GMA driver genes", 
                     top = 10)

# Check enriched pathways among the overlapping driver genes.
nrow(overlap_dr$data) # 0
