## This script performs DEA on gene expression data

## Load libraries ------------------------------------
library(stringr)
library(purrr)
library(dplyr)
library(magrittr)
library(TCGAbiolinks)
library(limma)

## Load data ------------------------------------

# Load filtered gene expression data subsetted to contain same patients as 
# methylation data
data_exp <- get(load("breast_basal/data/BRCA_basal_dataFilt_wrangled.rda"))


## Analyze data ------------------------------------

# Split barcodes of expression data into its components
split_barcode <- str_split(colnames(data_exp), 
                           "-")

# Get TSS IDs (batch factor) from barcode
TSS_ID <- split_barcode %>%
  map_chr(~ .[[2]])

# Access conditions (cancer vs normal) from barcode
condition_samples <- split_barcode %>%
  map_chr(~ .[[4]])

# Create design matrix containing condition types and TSS 
design_matrix <- model.matrix(~0+condition_samples+TSS_ID)

# Assign colnames
colnames(design_matrix)[1:2] <- c("Cancer", "Normal")

# Voom transform data
voom_data <- voom(counts = data_exp,
                  design = design_matrix)

# Make group contrasts
normal_cancer_contrast <- makeContrasts("Cancer-Normal", 
                                        levels = design_matrix)

### Do DEA with limma 

# Fit linear model to processed data incorporating design matrix
lm_fit <- lmFit(object = voom_data, 
                design = design_matrix)

# Compute estimated coefficients for contrasts
contr_fit <- contrasts.fit(fit = lm_fit, 
                           contrasts = normal_cancer_contrast)

# Fit empirical Bayes model for differential expression using treat
eBayes_fit <- eBayes(fit = contr_fit)

# Extract table of genes from fit
DEA_table <- topTable(fit = eBayes_fit, 
                      coef = 1, 
                      adjust.method = "fdr", 
                      number = nrow(voom_data))

# Get significant DEGs
index_up <- which(DEA_table$logFC >= 1 & DEA_table$adj.P.Val < 0.05)
index_down <- which(DEA_table$logFC <= -1 & DEA_table$adj.P.Val < 0.05)
index_up_down <- c(index_up, 
                   index_down)
DEA_table_extract <- DEA_table[index_up_down,]


## Save data ------------------------------------

# Save table of all DEGs 
save(DEA_table, 
     file = "breast_basal/data/BRCA_basal_full_DEA.rda")

# Save table of significant DEGs
save(DEA_table_extract,
     file = "breast_basal/data/BRCA_basal_sig_DEA.rda")
