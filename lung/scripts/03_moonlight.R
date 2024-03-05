## This script performs DEA on gene expression data

## Load libraries ------------------------------------
library(Moonlight2R)
library(TCGAbiolinks)
library(magrittr)

## Load data ------------------------------------

# Load filtered gene expression data subsetted to contain same patients as 
# methylation data
data_exp <- get(load("lung/data/LUAD_dataFilt_wrangled.rda"))

# Load table of DEGs
data_DEGs <- get(load("lung/data/LUAD_sig_DEA.rda"))

# Load data from Moonlight package
data(DiseaseList)
data(EAGenes)
data(tabGrowBlock)
data(knownDriverGenes)

## Wrangle data ------------------------------------

# Get tumor samples from gene expression matrix
tumor_samples <- TCGAquery_SampleTypes(barcode = colnames(data_exp), 
                                       typesample = "TP")

# Keep only tumor samples in gene expression matrix
data_exp_tumor <- data_exp[, tumor_samples]


## Analyze data ------------------------------------

## Apply Moonlight framework

# Do functional enrichment analysis
data_FEA <- FEA(DEGsmatrix = data_DEGs)

print("FEA is done")

# Do gene regulatory network analysis
set.seed(317)
data_GRN <- GRN(TFs = rownames(data_DEGs), 
                DEGsmatrix = data_DEGs, 
                DiffGenes = TRUE, 
                normCounts = data_exp_tumor, 
                kNearest = 3, 
                nGenesPerm = 2000, 
                nBoot = 400)

print("GRN is done")

# Do upstream regulatory analysis
data_URA <- URA(dataGRN = data_GRN, 
                DEGsmatrix = data_DEGs, 
                BPname = c("apoptosis", 
                           "proliferation of cells"),
                nCores = 4)

print("URA is done")

# Do pattern recognition analysis
data_PRA <- PRA(dataURA = data_URA, 
                BPname = c("apoptosis", "proliferation of cells"), 
                thres.role = 0)

print("PRA is done")


## Save data ------------------------------------

# Save FEA data
save(data_FEA, file = "lung/data/LUAD_dataFEA.rda")

# Save GRN data
save(data_GRN, file = "lung/data/LUAD_dataGRN.rda")

# Save URA data
save(data_URA, file = "lung/data/LUAD_dataURA.rda")

# Save PRA data
save(data_PRA, file = "lung/data/LUAD_dataPRA.rda")


