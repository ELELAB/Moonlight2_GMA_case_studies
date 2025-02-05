## This script compares results across cancer types

## Load libraries ------------------------------------
library(Moonlight2R)
library(viridis)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)

## Source functions  --------------------------------
source("99_functions.R")

## Load data ------------------------------------

# Load overlap percentages and sensitivities 
overlaps_breast_basal <- read_csv("breast_basal/results/07_gene_overlap_cgc.csv")
overlaps_lung <- read_csv("lung/results/07_gene_overlap_cgc.csv")
overlaps_thyroid <- read_csv("thyroid/results/07_gene_overlap_cgc.csv")

# Load driver genes 
drivers_breast_basal <- get(load("breast_basal/results/BRCA_basal_dataEDA.rda"))
drivers_lung <- get(load("lung/results/LUAD_dataEDA.rda"))
drivers_thyroid <- get(load("thyroid/results/THCA_dataEDA.rda"))

# Load PRA data
data_PRA_breast <- get(load("breast_basal/data/BRCA_basal_dataPRA.rda"))
data_PRA_lung <- get(load("lung/data/LUAD_dataPRA.rda"))
data_PRA_thyroid <- get(load("thyroid/data/THCA_dataPRA.rda"))

# Load results of Cox multivariate regression analysis
cox_multi_results_lung <- get(load("lung/results/09_LUAD_cox_multivariate.rda"))
cox_multi_results_thyroid <- get(load("thyroid/results/09_THCA_cox_multivariate.rda"))

# Load driver gene-drug target interactions for the three cancer types
basal_drug <- get(load("breast_basal/results/10_BRCA_basal_drug_gene_interactions.rda")) 
lung_drug <- get(load("lung/results/10_LUAD_drug_gene_interactions.rda"))
thyroid_drug <- get(load("thyroid/results/10_THCA_drug_gene_interactions.rda"))

# Load oncogenic mediator summaries for the three cancer types
om_basal <- get(load("breast_basal/results/Oncogenic_mediators_methylation_summary.rda"))
om_lung <- get(load("lung/results/Oncogenic_mediators_methylation_summary.rda"))
om_thyroid <- get(load("thyroid/results/Oncogenic_mediators_methylation_summary.rda"))


## Wrangle data ----------------------------------

###### FOR GENE QUANTITATIVE STATISTICS ###### 

# Combine overlap overviews from cancer types into one list
overlaps_list <- list(overlaps_breast_basal, overlaps_lung, overlaps_thyroid)
names(overlaps_list) <- c("breast_basal", "lung", "thyroid")

# Select sensitivity and precision from overlap overviews
# Add annotation with cancer type
# Convert to long format to prepare data for plotting
overlaps_cancer_types <- map(seq.int(length(overlaps_list)), function(x) {
  overlaps_list[[x]] %>% 
    dplyr::select(gene_set, precision, sensitivity) %>% 
    dplyr::mutate(cancer_type = names(overlaps_list)[[x]]) }) %>% 
  bind_rows() %>% 
  pivot_longer(data = ., cols = -c("gene_set", "cancer_type"),
               names_to = "estimate_type", values_to = "estimate_perc")


###### FOR DUAL ROLE INVESTIGATION ###### 

# Combine predicted TSGs and OCGs of cancer types into one list
names(drivers_breast_basal) <- c("TSG_breast_basal", "OCG_breast_basal")
names(drivers_lung) <- c("TSG_lung", "OCG_lung")
names(drivers_thyroid) <- c("TSG_thyroid", "OCG_thyroid")
drivers_list <- c(drivers_breast_basal, drivers_lung, drivers_thyroid)

# Prepare predicted TSGs and OCGs data of cancer types for circos plotting
data_PRA_cancers <- list(data_PRA_breast, data_PRA_lung, data_PRA_thyroid)
data_drivers <- list(drivers_breast_basal, drivers_lung, drivers_thyroid)
names(data_PRA_cancers) <- names(data_drivers) <- c("Basal-like BRCA",
                                                    "LUAD", "THCA")
circos_data <- map2(data_PRA_cancers, data_drivers, function(x, y) {
  x$TSG <- x$TSG[y$TSG]
  x$OCG <- x$OCG[y$OCG]
  list(listCandidates = list("TSG" = x$TSG, "OCG" = x$OCG))
})


###### FOR HAZARD RATIO FOREST PLOT ###### 

# Filter Cox multivariate results of driver genes
# to include only significant results
# and add column indicating the cancer type
# and combine results for the cancer types
cox_multi_results_sig_lung <- cox_multi_results_lung$driver_genes_OCGs %>% 
  dplyr::filter(variable == "exp_value" & `Pr(>|z|)` < 0.05) %>% 
  dplyr::mutate(cancer_type = "LUAD")
cox_multi_results_sig_thyroid <- cox_multi_results_thyroid$driver_genes_OCGs %>% 
  dplyr::filter(variable == "exp_value" & `Pr(>|z|)` < 0.05) %>% 
  dplyr::mutate(cancer_type = "THCA")
cox_multi_results_cancers <- cox_multi_results_sig_lung %>% 
  bind_rows(cox_multi_results_sig_thyroid) %>% 
  dplyr::rename("hazard_ratio" = "exp(coef)",
                "CI_95_lower" = "2.5 %",
                "CI_95_upper" = "97.5 %")


###### FOR DRUG TARGET PLOT ###### 

# Extract results showing number of drugs that each driver gene interacts with
# for all cancer types in one tibble
drug_count_cancers <- map2(list(basal_drug, lung_drug, thyroid_drug), 
                           c("basal", "lung", "thyroid"), function(x, y) {
                             x@byGene %>% 
                               as_tibble() %>% 
                               dplyr::filter(DistinctDrugCount > 0) %>% 
                               dplyr::select(DistinctDrugCount) %>% 
                               dplyr::mutate(cancer_type = y)
                           }) %>% 
  bind_rows() 


###### FOR OVERVIEW OF DIFFERENTIALLY METHYLATED CPGs ###### 

# Define methylation status types
met_status <- c("Hyper_methyl_num", "Hypo_methyl_num", "Dual_methyl_num")

# Combine oncogenic mediator summaries for three cancer types into one list
om_cancers <- list(om_basal, om_lung, om_thyroid)
names(om_cancers) <- c("basal", "lung", "thyroid")

# Extract number of hypo, hyper and dual differentially methylated CpGs
# associated with the oncogenic mediators for the three cancer types
# and combine into one table
om_cancers_methyl_num <- map(1:length(om_cancers), function(x) {
  om_cancers[[x]] %>% 
    dplyr::select(all_of(met_status)) %>% 
    summarise(across(everything(), ~ sum(.))) %>% 
    dplyr::mutate(cancer = names(om_cancers)[[x]])
}) %>% 
  bind_rows()



## Visualize data ----------------------------------

# Visualize precision and sensitivity between gene set and CGC (in percentage) across
# cancer types 
overlaps_cancer_types_plot <- ggplot(data = overlaps_cancer_types, 
                                     mapping = aes(x = gene_set, 
                                                   y = estimate_perc,
                                                   fill = estimate_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  facet_grid(~ cancer_type) + 
  theme_bw() +
  labs(x = "", y = "Estimate [%]",
       title = "Precision and sensitivity [%] of methods stratified by cancer type") +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks = element_line(size = 2),
        title = element_text(size = 14, color = "black"),
        panel.border = element_rect(size = 2),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14)) 

# Exclude EpiMix from sensitivity/precision visualization
overlaps_cancer_types_wo_EpiMix <- overlaps_cancer_types %>% 
  dplyr::filter(gene_set != "EpiMix") %>% 
  mutate(gene_set = ifelse(gene_set == "EDA", "GMA", gene_set))
overlaps_cancer_types_plot_wo_EpiMix <- ggplot(data = overlaps_cancer_types_wo_EpiMix, 
                                               mapping = aes(x = estimate_type, 
                                                             y = estimate_perc,
                                                             fill = gene_set)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  facet_grid(~ cancer_type) + 
  theme_bw() +
  labs(x = "", y = "Estimate [%]",
       title = "Precision and sensitivity [%] of methods stratified by cancer type") +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks = element_line(size = 2),
        title = element_text(size = 14, color = "black"),
        panel.border = element_rect(size = 2),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14)) 

# Visualize dual role driver genes across cancer types in UpSet plot
dual_genes <- upset_plot(file_name = "98_dual_genes_cancer_types.pdf", 
                         dir_output = ".", sets_list = drivers_list, 
                         names_sets = names(drivers_list), 
                         title_plot = "Dual role driver genes across three cancer types")

# Visualize comparison of predicted driver genes across cancer types in circos plot
plotCircos(listMoonlight = circos_data, intensityColOCG = 1, intensityColTSG = 1, 
           fontSize = 1.5, additionalFilename = "_cancer_types_98")

# Visualize hazard ratios of prognostic genes at multivariate level
# in forest plot stratified by cancer type
forest_cox_mutli <- cox_multi_results_cancers %>% 
  ggplot(data = ., mapping = aes(x = fct_reorder(gene, cancer_type, .desc = TRUE), 
                                 y = hazard_ratio,
                                 ymin = CI_95_lower, ymax = CI_95_upper, 
                                 col = cancer_type)) +
  geom_pointrange(lwd = 1.5, size = 1.5) +
  scale_color_manual(values = c("#440154FF", "#FDE725FF")) +
  coord_flip() +
  labs(x = "", title = "Hazard ratios of OCGs stratified by cancer type") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.ticks = element_line(size = 1),
        title = element_text(size = 14, color = "black"),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white", size = 1.5),
        strip.text = element_text(size = 18))

# Visualize number of driver gene-drug target interactions
# with number of driver genes on y axis and
# number of drugs that the genes interact with on the x axis
drug_count_cancers_plot <- drug_count_cancers %>% 
  ggplot(data = ., mapping = aes(x = DistinctDrugCount, fill = cancer_type)) +
  geom_histogram(position = "stack", binwidth = 1) +
  labs(x = "Number of drug interactions",
       y = "Number of driver genes",
       title = "Distribution of driver gene-drug target interactions stratified by cancer type") +
  scale_x_continuous(breaks = seq(0, 55, 2)) +
  scale_y_continuous(breaks = seq(0, 26, 2)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        title = element_text(size = 14, color = "black"),
        panel.border = element_rect(size = 0),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "cm")) 

# Visualize number of differentially methylated CpGs in the oncogenic 
# mediators in the three cancer types
om_cancers_methyl_num_plot <- om_cancers_methyl_num %>% 
  pivot_longer(cols = -c(cancer)) %>% 
  ggplot(., mapping = aes(x = name, y = value, fill = cancer)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 5) + 
  labs(x = "Methylation type", 
       y = "Number of differentially methylated CpGs in oncogenic mediators", 
       fill = "Cancer type") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        panel.border = element_rect(linewidth = 1, colour = "black"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14)) 


## Save data ------------------------------------

# Save gene overlap plots
ggsave(filename = "98_prec_sens_cancer_types_epimix_moonlight.pdf",
       plot = overlaps_cancer_types_plot,
       path = ".",
       width = 10,
       height = 5)
ggsave(filename = "98_prec_sens_cancer_types_moonlight.pdf",
       plot = overlaps_cancer_types_plot_wo_EpiMix,
       path = ".",
       width = 10,
       height = 5)

# Save forest plot of hazard ratios of driver genes 
ggsave(filename = "98_forest_plot_cox_multi_cancers.pdf",
       plot = forest_cox_mutli,
       path = ".",
       width = 8,
       height = 6)

# Save distribution plot of driver gene-drug target interactions
# across all three cancer types
ggsave(filename = "98_driver_gene_drugs_cancers.pdf",
       plot = drug_count_cancers_plot,
       path = ".",
       width = 10,
       height = 4)

# Save barplot showing number of differentially methylated CpGs in the oncogenic 
# mediators in the three cancer types
ggsave(filename = "98_diff_met_cpgs_oncogenic_mediators.pdf",
       plot = om_cancers_methyl_num_plot,
       path = ".",
       width = 10,
       height = 4)
