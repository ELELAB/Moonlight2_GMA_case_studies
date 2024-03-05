## This script visualizes results from Moonlight and EpiMix

## Load libraries ------------------------------------
library(Moonlight2R)
library(RColorBrewer)
library(VennDiagram)
library(tidyverse)


## Load data ------------------------------------

# Load FEA data
data_FEA <- get(load("breast_basal/data/BRCA_basal_dataFEA.rda"))

# Load DEA + methylation overview table
DEG_Methylation_Annotations <- get(load("breast_basal/results/DEG_Methylation_Annotations.rda"))

# Load oncogenic mediators + number of methylated CpGs overview table
Oncogenic_mediators_methylation_summary <- get(load("breast_basal/results/Oncogenic_mediators_methylation_summary.rda"))

# Load URA data
data_URA <- get(load("breast_basal/data/BRCA_basal_dataURA.rda"))

# Load PRA data (Moonlight's primary layer)
data_PRA <- get(load("breast_basal/data/BRCA_basal_dataPRA.rda"))

# Load functional gene-CpG pairs (EpiMix)
functional_pairs <- read_csv("breast_basal/results/FunctionalPairs_Regular.csv")

# Load GMA data (list of driver genes)
data_EDA <- get(load("breast_basal/results/BRCA_basal_dataEDA.rda"))


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


## Visualize data ----------------------------------

## Visualize FEA

# Find significantly enriched BPs from FEA
# Biological processes are significantly enriched if 
# absolute value of Moonlight Z-score >= 1 and FDR <= 0.01
data_FEA_sig <- data_FEA %>% 
  dplyr::filter(abs(Moonlight.Z.score) >= 1, 
                FDR <= 0.01)

# Visualize FEA results
plotFEA(dataFEA = data_FEA_sig, 
        topBP = nrow(data_FEA_sig), 
        additionalFilename = "breast_basal/results/05_BRCA_basal_FEAplot.pdf",
        height = 6, 
        width = 10)

## Visualize number of hypo/hyper/dual methylated CpG sites in oncogenic mediators

# plotGMA function with complete mode
plotGMA(DEG_Methylation_Annotations = DEG_Methylation_Annotations, 
        Oncogenic_mediators_methylation_summary = Oncogenic_mediators_methylation_summary, 
        type = "complete", 
        additionalFilename = "breast_basal/results/05_BRCA_basal_")

## Visualize effect of driver genes on biological processes

# plotMoonlightMet function 
plotMoonlightMet(DEG_Methylation_Annotations = DEG_Methylation_Annotations, 
                 Oncogenic_mediators_methylation_summary = Oncogenic_mediators_methylation_summary, 
                 dataURA = data_URA, 
                 genes = driver_genes, 
                 additionalFilename = "breast_basal/results/05_BRCA_basal_")

## Visualize predicted genes from Moonlight's primary layer, EpiMix and GMA

# Define colors to use in Venn diagram
venn.cols <- brewer.pal(3, "Pastel2")

# Create Venn diagram 
venn.diagram(gene_lists, filename = "breast_basal/results/05_BRCA_basal_venn_gene_lists.png", 
             category.names = c("Moonlight's \n primary \n layer", "EpiMix", "GMA"), output = TRUE, 
             main = "Venn diagram of gene lists", 
             main.cex = 1.5, main.fontface = "bold", 
             main.fontfamily = "sans", lwd = 3, 
             lty = "blank", fill = venn.cols, 
             cex = 2, fontface = "bold", 
             fontfamily = "sans", cat.cex = 1.5, 
             cat.fontface = "bold", cat.fontfamily = "sans",
             cat.dist = c(0.00005, 0.00005, 0.00005),
             cat.pos = c(320, -320, -180))



