## This script performs enrichment analysis on the DMA predicted driver genes and checks for overlaps between these drivers and the GMA predicted drivers

## Load libraries -----------
library(Moonlight2R)
library(enrichR)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(VennDiagram)

## Source functions script  --------------------------
source("99_functions.R")

## Load data ---------------
load("lung/results/LUAD_dataEDA.rda")
load("lung/results/LUAD_dataDMA.rda")

# Get driver genes
GMA_driver_genes <- Reduce(c, data_EDA)
DMA_driver_genes <- Reduce(c, data_DMA)

# Wrangle the two gene lists into one list
gene_lists <- list(GMA_driver_genes,DMA_driver_genes)
names(gene_lists) <- c('GMA_driver','DMA_driver')

# Separate GMA TSGs and OCGs into lists
GMA_TSG <- data_EDA$TSG
GMA_OCG <- data_EDA$OCG
DMA_TSG <- data_DMA$TSG
DMA_OCG <- data_DMA$OCG

## Analyze data ---------

## Venn diagrams --------

# Define colors to use in Venn diagram
venn.cols <- c("#40E0D0", "#C8A2C8")

# Make diagrams for comparison of TSGs, OCGs and all genes between the two functions
gene_list <- list(GMA_TSG, DMA_TSG)
venn.diagram(gene_list, filename = "lung/results/14_LUAD_venn_GMA_DMA_TSG.png", 
             category.names = c("GMA",  "DMA"), output = TRUE, 
             main = "Count of TSGs predicted by GMA and DMA respectively", 
             main.cex = 1.2, main.fontface = "bold", 
             main.fontfamily = "sans", lwd = 2, 
             lty = "blank", fill = venn.cols, 
             cex = 2, fontface = "bold", 
             fontfamily = "sans", cat.cex = 1.5, 
             cat.fontface = "bold", cat.fontfamily = "sans",
             cat.dist = c(0.00005, 0.00005),
             cat.pos = c( 200, -200))

gene_list <- list(GMA_OCG, DMA_OCG)
venn.diagram(gene_list, filename = "lung/results/14_LUAD_venn_GMA_DMA_OCG.png", 
             category.names = c("GMA",  "DMA"), output = TRUE, 
             main = "Count of OCGs predicted by GMA and DMA respectively", 
             main.cex = 1.2, main.fontface = "bold", 
             main.fontfamily = "sans", lwd = 2, 
             lty = "blank", fill = venn.cols, 
             cex = 2, fontface = "bold", 
             fontfamily = "sans", cat.cex = 1.5, 
             cat.fontface = "bold", cat.fontfamily = "sans",
             cat.dist = c(0.00005, 0.00005),
             cat.pos = c( 200, -200), ext.text = FALSE)

gene_list <- list(GMA_driver_genes, DMA_driver_genes)
venn.diagram(gene_list, filename = "lung/results/14_LUAD_venn_GMA_DMA_all.png", 
             category.names = c("GMA",  "DMA"), output = TRUE, 
             main = "Count of driver genes predicted by GMA and DMA respectively", 
             main.cex = 1.1, main.fontface = "bold", 
             main.fontfamily = "sans", lwd = 2, 
             lty = "blank", fill = venn.cols, 
             cex = 2, fontface = "bold", 
             fontfamily = "sans", cat.cex = 1.5, 
             cat.fontface = "bold", cat.fontfamily = "sans",
             cat.dist = c(0.00005, 0.00005),
             cat.pos = c( 200, -200))


## Perform enrichment analysis -----------

# Load databases for enrichment analysis
dbs <- c("MSigDB_Hallmark_2020")

# Perform enrichment analysis on each of the gene lists
enrich_analysis <- map(gene_lists, function(x) { 
  enrichr(genes = x, databases = dbs)
})

## Visualize data -----------

# Visualize enrichment analysis results of MSigDB Hallmark terms of GMA driver genes
gma_dr <- goplot(data = enrich_analysis$GMA_driver$MSigDB_Hallmark_2020, 
                title = "MSigDB Hallmark 2020 - GMA driver genes", 
                top = 10)

dma_dr <- goplot(data = enrich_analysis$DMA_driver$MSigDB_Hallmark_2020, 
                  title = "MSigDB Hallmark 2020 - DMA driver genes", 
                  top = 10)

# Visualize enrichment analysis of all gene lists side by side
all_dr <- gma_dr + dma_dr

# Save side by side enrichment plot of MSigDB hallmark terms
ggsave(filename = "14_LUAD_GMA_vs_DMA.pdf",
       plot = all_dr,
       path = "lung/results/",
       width = 12,
       height = 8)
