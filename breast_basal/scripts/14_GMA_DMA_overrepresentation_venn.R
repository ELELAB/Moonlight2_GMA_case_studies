## This script performs enrichment analysis on the DMA predicted driver genes and checks for overlaps between these drivers and the GMA predicted drivers

## Load libraries ------------------------------------
library(Moonlight2R)
library(enrichR)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(VennDiagram)

## Source functions script  --------------------------
source("99_functions.R")

## Load the GMA and DMA data  --------------------------
data_EDA <- get(load("breast_basal/results/BRCA_basal_dataEDA.rda"))
data_DMA <- get(load("breast_basal/data/BRCA_basal_dataDMA.rda"))

## Get driver genes  --------------------------
GMA_driver_genes <- Reduce(c, data_EDA)
DMA_driver_genes <- Reduce(c, data_DMA)

## Wrangle the two gene lists into one list
gene_lists <- list(GMA_driver_genes,DMA_driver_genes)
names(gene_lists) <- c('GMA_driver','DMA_driver')


## Analyze data ----------------------------------

## Venn Diagram
GMA_TSG <- data_EDA$TSG
GMA_OCG <- data_EDA$OCG

DMA_TSG <- DMA_Basal$TSG
DMA_OCG <- DMA_Basal$OCG

# Define colors to use in Venn diagram
venn.cols <- c("#5bacf3", "#e4f62f")

# Create Venn diagram 
gene_list <- list(GMA_TSG,DMA_TSG)
venn.diagram(gene_list, filename = "breast_basal/results/14_BRCA_basal_venn_GMA_DMA_TSG.png", 
             category.names = c("GMA",  "DMA"), output = TRUE, 
             main = "Count of TSGs predicted by GMA and DMA respectively ", 
             main.cex = 1.2, main.fontface = "bold", 
             main.fontfamily = "sans", lwd = 2, 
             lty = "blank", fill = venn.cols, 
             cex = 2, fontface = "bold", 
             fontfamily = "sans", cat.cex = 1.5, 
             cat.fontface = "bold", cat.fontfamily = "sans",
             cat.dist = c(0.00005, 0.00005),
             cat.pos = c( 200, -200))

# Create Venn diagram 
gene_list <- list(GMA_OCG,DMA_OCG)
venn.diagram(gene_list, filename = "breast_basal/results/14_BRCA_basal_venn_GMA_DMA_OCG.png", 
             category.names = c("GMA",  "DMA"), output = TRUE, 
             main = "Count of OCGs predicted by GMA and DMA respectively ", 
             main.cex = 1.2, main.fontface = "bold", 
             main.fontfamily = "sans", lwd = 2, 
             lty = "blank", fill = venn.cols, 
             cex = 2, fontface = "bold", 
             fontfamily = "sans", cat.cex = 1.5, 
             cat.fontface = "bold", cat.fontfamily = "sans",
             cat.dist = c(0.00005, 0.00005),
             cat.pos = c( 320, 320))

# Create Venn diagram 
gene_list <- list(GMA_driver_genes,DMA_driver_genes)
venn.diagram(gene_list, filename = "breast_basal/results/14_BRCA_basal_venn_GMA_DMA_all.png", 
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

## Perform enrichment analysis
enrich_analysis <- map(gene_lists, function(x) { 
  enrichr(genes = x, databases = "MSigDB_Hallmark_2020")
})

## Visualize data ----------------------------------

# Visualize enrichment analysis results of MSigDB Hallmark terms of driver genes
GMA_Hm <-goplot(data = enrich_analysis$GMA_driver$MSigDB_Hallmark_2020,
                title = "MSigDB Hallmark 2020 - GMA driver genes",
                top = 10)

DMA_Hm <-goplot(data = enrich_analysis$DMA_driver$MSigDB_Hallmark_2020,
                title = "MSigDB Hallmark 2020 - DMA driver genes",
                top = 10)

# Visualize enrichment analysis of all gene lists side by side
all_Hm <- GMA_Hm + DMA_Hm

## Save data ------------------------------------
ggsave(filename = "14_comparison_GMA_DMA.pdf",
       plot = all_Hm,
       path = "breast_basal/results/",
       width = 12,
       height = 8)
