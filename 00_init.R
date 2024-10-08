### This script runs all analyses associated with the discovery of methylation-
### driven driver genes in cancer (sub)types using the Moonlight2R framework.
### The data that is analyzed is from The Cancer Genome Atlas.


# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Download data from OSF --------------------------------------------------
#library(osfr)
#node <- osf_retrieve_node("j4n8q")
#files <- osf_ls_files(node, n_max = Inf)
#osf_download(files, 
#             recurse = TRUE, 
#             conflicts = "skip")

# Restore packages from lockfile -----------------------------------------
library(renv)
renv::settings$external.libraries(c("methyl_case/lib/R/library/"))
renv::restore()
renv::activate()
detach("package:renv", unload = TRUE)

# Set up for Cairo rendering, for headless machines ----------------------
options(bitmapType = "cairo"


## Run scripts ------------------------------------------------------------

## Lung adenocarcinoma -----------------------------------------------------

# Get and wrangle expression and methylation data
source("lung/scripts/01_data_wrangling.R")

# Perform differential expression analysis between normal and cancer tissue
source("lung/scripts/02_DEA.R")

# Run Moonlight's primary layer
source("lung/scripts/03_moonlight.R")

# Run Moonlight's secondary methylation layer
source("lung/scripts/04_epimix_moonlight.R")

# Visualize results from Moonlight
source("lung/scripts/05_epimix_moonlight_visualization.R")

# Perform enrichment analysis of predicted driver genes
source("lung/scripts/06_enrich_analysis.R")

# Compare predicted driver genes with the Cancer Gene Census
source("lung/scripts/07_cancer_gene_census.R")

# Perform survival analysis of predicted driver genes
source("lung/scripts/09_survival_analysis.R")

# Investigate drug targets among predicted driver genes
source("lung/scripts/10_drug_gene_interactions.R")

# Subset mutation samples from methylation dataset
source("lung/scripts/11_mutations_subset.R")

# Predict DMA driver genes from the mutation subset
source("lung/scripts/12_DMA_lung.R")

# Check for samples with driver mutations but no methylation
source("lung/scripts/13_filtering_nonGMA.R")

# Perform enrichment analysis on DMA predicted driver genes and look for overlaps with GMA
source("lung/scripts/14_GMA_DMA_overrepresentation_venn.R")

# Perform enrichment analysis on overlapping driver genes between GMA and DMA
source("lung/scripts/15_GMA_DMA_overlap_enrichment.R")

## Compare results from the cancer (sub)types -----------------------------

# Compare results from the cancer (sub)types
# source("98_compare_cancer_types.R")
