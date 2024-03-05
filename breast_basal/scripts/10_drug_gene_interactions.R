## This script finds investigates drug-gene interactions
## among the predicted driver genes

## Load libraries ------------------------------------
library(tidyverse)
library(rDGIdb)


## Load data ------------------------------------

# Load driver genes (Moonlight + EpiMix i.e. Moonlight's secondary layer)
data_EDA <- get(load("breast_basal/results/BRCA_basal_dataEDA.rda"))


## Wrangle data ----------------------------------

# Get driver genes
driver_genes <- Reduce(c, data_EDA)


## Analyze data ----------------------------------

# Query DGIdb for driver gene-drug interactions
# using only cancer-specific source databases
drug_gene_interact <- queryDGIdb(driver_genes, 
                                 sourceDatabases = c("CGI", "CIViC", "COSMIC",
                                                     "CancerCommons", 
                                                     "ClearityFoundationBiomarkers",
                                                     "ClearityFoundationClinicalTrial",
                                                     "DoCM", "JAX-CKB", "MyCancerGenome",
                                                     "MyCancerGenomeClinicalTrial",
                                                     "NCI", "OncoKB", "TALC"))

# Extract table showing number of drugs each driver gene interacts with
# Add annotation if driver gene is TSG or OCG
drug_gene_count <- drug_gene_interact@byGene %>% 
  as_tibble() %>% 
  dplyr::mutate(driver_type = case_when(Gene %in% data_EDA$TSG ~ "TSG",
                                        Gene %in% data_EDA$OCG ~ "OCG"))

# Extract table showing which drugs interact with which driver genes
# Add annotation if driver gene is TSG or OCG
# Add annotation if gene and drug interact (yes = interaction)
drug_gene_names <- drug_gene_interact@detailedResults %>% 
  as_tibble() %>% 
  dplyr::mutate(driver_type = case_when(Gene %in% data_EDA$TSG ~ "TSG",
                                        Gene %in% data_EDA$OCG ~ "OCG"),
                Interaction = "Yes")

# Create tibble showing all possible combinations of driver gene-drug
# interactions (only for those driver genes that interact with
# at least 1 drug) 
all_genes <- drug_gene_names$Gene %>% unique()
all_drugs <- drug_gene_names$Drug %>% unique()
all_combinations <- expand.grid(Gene = all_genes, Drug = all_drugs) %>% 
  as_tibble()

# Join driver gene-drug tibble with all combinations tibble
# where yes indicates driver gene and drug interact and 
# no indicates driver gene and drug do not interact
drug_gene_names_full <- full_join(x = all_combinations, y = drug_gene_names, 
                                  by = c("Gene", "Drug")) %>% 
  dplyr::select(c("Gene", "Drug", "Interaction", "InteractionType")) %>% 
  dplyr::mutate(driver_type = case_when(Gene %in% data_EDA$TSG ~ "TSG",
                                        Gene %in% data_EDA$OCG ~ "OCG"),
                Interaction = case_when(Interaction == "Yes" ~ "Yes",
                                        is.na(Interaction) ~ "No"),
                Interaction = as.factor(Interaction))

# Create same tibble as above but containing only those
# driver gene-drug target combinations where the interaction
# type is known
drug_gene_names_full_type <- drug_gene_names_full %>% 
  dplyr::filter(!is.na(InteractionType),
                !InteractionType == "")
all_genes_type <- drug_gene_names_full_type$Gene %>% unique()
all_drugs_type <- drug_gene_names_full_type$Drug %>% unique()
all_combinations_type <- expand.grid(Gene = all_genes_type, Drug = all_drugs_type) %>% 
  as_tibble()
drug_gene_names_full_type <- full_join(x = all_combinations_type, y = drug_gene_names_full_type,
                                       by = c("Gene", "Drug")) %>% 
  dplyr::mutate(driver_type = case_when(Gene %in% data_EDA$TSG ~ "TSG",
                                        Gene %in% data_EDA$OCG ~ "OCG"),
                Interaction = case_when(Interaction == "Yes" ~ "Yes",
                                        is.na(Interaction) ~ "No"),
                Interaction = as.factor(Interaction),
                InteractionType = case_when(is.na(InteractionType) ~ "No",
                                            TRUE ~ InteractionType))



## Visualize data ----------------------------------

# Visualize number of drugs each driver gene interacts with
# Show only results for those driver genes that interact with
# at least 1 drug
drug_gene_count_data_plot <- drug_gene_count %>% 
  dplyr::filter(DistinctDrugCount > 0) 
cols_groups <- case_when(drug_gene_count_data_plot$driver_type %>% 
                           unique %>% 
                           length > 1 ~ c("#21908CFF", "#440154FF"),
                         TRUE ~ "#21908CFF")
drug_gene_count_plot <- ggplot(data = drug_gene_count_data_plot, 
                               mapping = aes(x = DistinctDrugCount, 
                                             y = reorder(Gene, DistinctDrugCount),
                                             fill = driver_type)) +
  geom_col() +
  geom_text(mapping = aes(label = DistinctDrugCount), 
            position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2,
            size = 4, col = "white") +
  labs(x = "Number of drug interactions",
       y = "", title = "Number of driver gene-drug interactions") +
  scale_fill_manual(values = cols_groups) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18,
                                   color = "black"),
        axis.title.x = element_text(size = 18,
                                    color = "black"),
        axis.title.y = element_text(size = 18,
                                    color = "black"),
        axis.text.y = element_text(size = 18,
                                   color = "black"),
        title = element_text(size = 18,
                             color = "black"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.3, "cm"),
        panel.border = element_blank())

# Visualize driver gene-drug interactions in heatmap
# only for those driver genes that interact with at least 1 drug
drug_gene_heatmap <- drug_gene_names_full %>% 
  dplyr::mutate(Drug = case_when(str_length(Drug) > 30 ~ str_c(str_sub(Drug, start = 1, end = 14),
                                                               "-\n", str_sub(Drug, start = 14 + 1)),
                                 TRUE ~ Drug)) %>% 
  ggplot(., mapping = aes(x = Drug, y = Gene, 
                          fill = Interaction)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(.~ driver_type, scales = "free_y", dir = "v", strip.position = "right") +
  scale_x_discrete(position = "top", guide = guide_axis(angle = 90)) +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = c("grey", "#21908CFF")) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(line = element_blank(),
        axis.text.x = element_text(size = 12,
                                   color = "black"),
        axis.text.y = element_text(size = 12,
                                   color = "black"),
        title = element_text(size = 18,
                             color = "black"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "cm"))

# Visualize driver gene-drug interactions in heatmap
# only for those interactions where the interaction
# type is known 
drug_gene_type_heatmap <- drug_gene_names_full_type %>% 
  dplyr::mutate(Drug = case_when(str_length(Drug) > 30 ~ str_c(str_sub(Drug, start = 1, end = 14),
                                                               "-\n", str_sub(Drug, start = 14 + 1)),
                                 TRUE ~ Drug)) %>% 
  ggplot(., mapping = aes(x = Drug, y = Gene, 
                          fill = InteractionType)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(.~ driver_type, scales = "free_y", dir = "v", strip.position = "right") +
  scale_x_discrete(position = "top", guide = guide_axis(angle = 90)) +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = c("#000004FF", "#21908CFF", "grey")) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(line = element_blank(),
        axis.text.x = element_text(size = 12,
                                   color = "black"),
        axis.text.y = element_text(size = 12,
                                   color = "black"),
        title = element_text(size = 18,
                             color = "black"),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"))




## Save data ------------------------------------

# Save results of drug-gene interaction query
save(drug_gene_interact, file = "breast_basal/results/10_BRCA_basal_drug_gene_interactions.rda")

# Save barplot showing number of drug-gene interactions of driver genes
ggsave(filename = "10_BRCA_basal_drug_gene_interactions_num.pdf",
       plot = drug_gene_count_plot,
       path = "breast_basal/results/",
       width = 8,
       height = 4)

# Save heatmap showing drug-gene interactions of driver genes
ggsave(filename = "10_BRCA_basal_drug_gene_interactions_heatmap.pdf",
       plot = drug_gene_heatmap,
       path = "breast_basal/results/",
       width = 12,
       height = 4)

# Save heatmap showing drug-gene interactions of driver genes
ggsave(filename = "10_BRCA_basal_drug_gene_interactions_type_heatmap.pdf",
       plot = drug_gene_type_heatmap,
       path = "breast_basal/results/",
       width = 10,
       height = 2.6)







