## This script contains various functions used in the case studies

### ---------------------------------------------------------------------------------------------------------

# This function, goplot, visualizes results of an enrichment analysis.
# This function takes three arguments: 
# data:           table containing results from an enrichment analysis.
# title:          string specifying the title of the plot.
# top:            integer specifying number of terms to plot on the y-axis.
# The function returns a ggplot visualizing the enrichment analysis. 
goplot <- function(data, 
                   title, 
                   top) {
  myplot <- data %>% 
    separate(col = Overlap, 
             into = c("Count",
                      "Total"), 
             sep ="/", 
             remove = FALSE) %>% 
    mutate(Count = as.numeric(Count),
           Gene.Ratio = Count/as.numeric(Total)) %>% 
    filter(Adjusted.P.value < 0.05) %>%
    slice_min(order_by = Adjusted.P.value,
              n = top,
              with_ties = FALSE) %>% 
    ggplot(data = .,
           mapping = aes(x = Gene.Ratio,
                         y = reorder(Term, 
                                     Gene.Ratio), 
                         color = Adjusted.P.value, 
                         size = Count)) +
    geom_point() +
    scale_y_discrete(labels = function(Term) str_wrap(Term, width = 20)) +
    #scale_x_continuous(breaks = seq(from = 0,
    #                                to = 1,
    #                                by = 0.1)) +
    scale_color_gradient(low = "red", 
                         high = "blue") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12,
                                     color = "black"),
          axis.title.x = element_text(size = 12,
                                      color = "black"),
          axis.text.y = element_text(size = 12,
                                     color = "black"),
          axis.ticks = element_line(size = 2),
          title = element_text(size = 12,
                               color = "black"),
          panel.border = element_rect(size = 2),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.3, "cm")) +
    labs(x = "Gene ratio",
         y = "",
         title = title)
  return(myplot)
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, upset_plot, creates an UpSet plot of intersections between
### different sets. This function takes five inputs:
### file_name:      File name of UpSet plot to be saved as a pdf
### dir_output:     Directory where UpSet plot should be saved
### sets_list:      A list containing those sets that are to be compared 
### names_sets:     A character vector containing names of each element in the list, sets_list
### title_plot:     The title to be included in the UpSet plot
### The function returns an UpSet plot. 
upset_plot <- function(file_name,
                       dir_output,
                       sets_list, 
                       names_sets, 
                       title_plot) {
  
  # Save UpSet plot as pdf
  pdf(file = paste0(dir_output, "/", file_name),
      width = 18,
      height = 8)
  
  # Add names to elements in list of sets
  names(sets_list) <- names_sets
  
  # Make combination matrix of sets to be used in UpSet plot
  comb_mat_sets <- make_comb_mat(sets_list, mode = "intersect")
  
  # Generate color palette for plot using viridis package
  n_col <- max(comb_degree(comb_mat_sets))
  palette_col <- viridis_pal(option = "viridis")(n_col)
  
  # Create UpSet plot
  upset_p <- UpSet(comb_mat_sets, 
                   set_order = names(sets_list),
                   pt_size = unit(5, "mm"), 
                   lwd = 3, 
                   height = unit(4, "cm"),
                   comb_col = palette_col[comb_degree(comb_mat_sets)],
                   top_annotation = upset_top_annotation(comb_mat_sets, 
                                                         height = unit(12, "cm"),
                                                         ylim = c(0, max(comb_size(comb_mat_sets))),
                                                         bar_width = 0.7, 
                                                         axis_param = list(side = "left", 
                                                                           at = seq(from = 0,
                                                                                    to = max(comb_size(comb_mat_sets)),
                                                                                    by = 500)),
                                                         annotation_name_side = "left", 
                                                         annotation_name_gp = gpar(cex = 1), 
                                                         annotation_name_offset = unit(1.5, "cm")),
                   right_annotation = upset_right_annotation(comb_mat_sets, 
                                                             width = unit(3, "cm"), 
                                                             gp = gpar(fill = "darkseagreen"),
                                                             axis_param = list(at = seq(from = 0,
                                                                                        to = max(set_size(comb_mat_sets)),
                                                                                        by = 2000)), 
                                                             annotation_name_offset = unit(1.5, "cm")),
                   row_names_gp = gpar(fontsize = 12))
  
  # Add number of elements in each set on top of bars in plot
  draw_upset <- draw(upset_p)
  col_ord <- column_order(draw_upset)
  c_s <- comb_size(comb_mat_sets)
  decorate_annotation("intersection_size", {
    grid.text(c_s[col_ord], 
              x = seq(c_s), 
              y = unit(c_s[col_ord], "native") + 
                unit(2, "pt"), 
              gp = gpar(fontsize = 12, fontface = "bold"),
              just = "bottom",
              default.units = "native")
  })
  
  # Add title to plot
  grid.text(label = title_plot, 
            x = unit(20, "cm"), 
            y = unit(18, "cm"), 
            gp = gpar(fontsize = 18),
            just = "centre")
  
  # End with dev.off() to save plot
  dev.off()
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, cox_reg_ph, tests the proportional hazards assumption of
### Cox regression. 
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: stage containing the stage annotation of 
###                 each patient, time_years with the survival time, 
###                 vital_status_binary with the vital status encoded as 0/1,
###                 age_years with the age of the patients in years, and 
###                 gender_binary indicating the gender of the patients
###                 encoded as 0/1. 
### The function returns a tibble of the results of the proportional hazards
### assumption test. 
cox_reg_ph <- function(exp_data, gene, clinical_data) {
  
  ## Fit a Cox regression model for each gene
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "submitter_id") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  # and remove samples annotated with Stage X og NA in stage
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "submitter_id") %>% 
    dplyr::filter(!is.na(stage),
                  !stage == "Stage X",
                  !is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Fit a univariate Cox regression model with expression values of respective gene
  # as explanatory variable
  cox_fit_uni <- coxph(Surv(time_years, vital_status_binary)~exp_value, data = clinical_exp_gene)
  
  # Fit a multivariate Cox regression model with expression values, tumor stage,
  # age, and gender as explanatory variables 
  # If gender has only one unique value, do not include it in multivariate analysis
  if (length(unique(clinical_exp_gene$gender_binary)) == 1) {
    cox_fit_multi <- coxph(Surv(time_years, vital_status_binary)~exp_value+stage+age_years, data = clinical_exp_gene)  
  } else {
    cox_fit_multi <- coxph(Surv(time_years, vital_status_binary)~exp_value+stage+age_years+gender_binary, data = clinical_exp_gene)
  }
  
  # Check the proportional hazards assumption of the univariate analysis
  ph_check_uni <- cox.zph(cox_fit_uni)
  ph_check_uni <- ph_check_uni$table %>% 
    as_tibble(rownames = "test") %>% 
    dplyr::mutate("test_type" = "univariate")
  
  # Check the proportional hazards assumption of the multivariate analysis
  ph_check_multi <- cox.zph(cox_fit_multi)
  ph_check_multi <- ph_check_multi$table %>% 
    as_tibble(rownames = "test") %>% 
    dplyr::mutate("test_type" = "multivariate")
  
  # Bind results from univariate and multivariate analyses together
  ph_check_results <- ph_check_uni %>% bind_rows(ph_check_multi) %>% 
    dplyr::mutate("gene" = gene)
  
  
  return(ph_check_results)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, cox_reg_univariate, performs univariate Cox regression.
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: time_years with the survival time and 
###                 vital_status_binary with the vital status encoded as 0/1
### The function returns a tibble of the results of univariate 
### Cox regression analysis. 
cox_reg_univariate <- function(exp_data, gene, clinical_data) {
  
  ## Fit a univariate Cox regression model for a gene
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    BiocGenerics::as.data.frame() %>% 
    rownames_to_column(var = "submitter_id") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  # and remove samples annotated with Stage X og NA in stage
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "submitter_id") %>% 
    dplyr::filter(!is.na(stage),
                  !stage == "Stage X",
                  !is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Fit a univariate Cox regression model with expression values of respective gene
  # as explanatory variable
  cox_fit_uni <- coxph(Surv(time_years, vital_status_binary)~exp_value, data = clinical_exp_gene)
  
  # Extract the coefficients and p value
  sum_cox_uni <- summary(cox_fit_uni)$coefficients
  sum_cox_uni <- sum_cox_uni %>% as_tibble(rownames = "variable") %>% 
    dplyr::mutate("cox_model_type" = "univariate",
                  "gene" = gene)
  sum_cox_uni_95CI <- exp(confint(cox_fit_uni)) %>% 
    as_tibble(rownames = "variable") %>% 
    full_join(sum_cox_uni, by = "variable") %>% 
    dplyr::relocate(., c(`2.5 %`, `97.5 %`), .after = `Pr(>|z|)`)
  #cox_p_value_uni <- sum_cox_uni$`Pr(>|z|)`
  
  return(sum_cox_uni_95CI)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, cox_reg_multivariate, performs multivariate Cox regression.
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: stage containing the stage annotation of 
###                 each patient, time_years with the survival time, 
###                 vital_status_binary with the vital status encoded as 0/1,
###                 age_years with the age of the patients in years, and 
###                 gender_binary indicating the gender of the patients 
###                 encoded as 0/1. 
### The function returns a tibble of the results of multivariate 
### Cox regression analysis. 
cox_reg_multivariate <- function(exp_data, gene, clinical_data) {
  
  ## Fit a multivariate Cox regression model for a gene
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    BiocGenerics::as.data.frame() %>% 
    rownames_to_column(var = "submitter_id") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  # and remove samples annotated with Stage X og NA in stage
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "submitter_id") %>% 
    dplyr::filter(!is.na(stage),
                  !stage == "Stage X",
                  !is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Fit a multivariate Cox regression model with expression values, tumor stage,
  # age, and gender as explanatory variables 
  # If gender has only one unique value, do not include it in multivariate analysis
  if (length(unique(clinical_exp_gene$gender_binary)) == 1) {
    cox_fit_multi <- coxph(Surv(time_years, vital_status_binary)~exp_value+stage+age_years, data = clinical_exp_gene)  
  } else {
    cox_fit_multi <- coxph(Surv(time_years, vital_status_binary)~exp_value+stage+age_years+gender_binary, data = clinical_exp_gene)
  }
  
  # Extract the coefficients, hazard ratio, p value and 95% CI of the hazard ratio
  sum_cox_multi <- summary(cox_fit_multi)$coefficients
  sum_cox_multi <- sum_cox_multi %>% as_tibble(rownames = "variable") %>% 
    dplyr::mutate("cox_model_type" = "multivariate",
                  "gene" = gene)
  sum_cox_multi_95CI <- exp(confint(cox_fit_multi)) %>% 
    as_tibble(rownames = "variable") %>% 
    full_join(sum_cox_multi, by = "variable") %>% 
    dplyr::relocate(., c(`2.5 %`, `97.5 %`), .after = `Pr(>|z|)`)
  #cox_p_value_multi <- sum_cox_multi$`Pr(>|z|)`
  
  return(sum_cox_multi_95CI)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, km_survival, performs Kaplan-Meier survival analysis. 
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: stage containing the stage annotation of 
###                 each patient, time_years with the survival time, 
###                 vital_status_binary with the vital status encoded as 0/1,
###                 age_years with the age of the patients in years, and 
###                 gender indicating the gender of the patients. 
### The function returns the fitted KM survival object. 
km_survival <- function(exp_data, gene, clinical_data) {
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    BiocGenerics::as.data.frame() %>% 
    rownames_to_column(var = "submitter_id") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  # and remove samples annotated with Stage X og NA in stage
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "submitter_id") %>% 
    dplyr::filter(!is.na(stage),
                  !stage == "Stage X",
                  !is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Get median expression level of gene
  median_exp <- median(clinical_exp_gene$exp_value)
  
  # Divide samples into two expression groups:
  # one group with expression values > median expression
  # one group with expression values < median expression
  clinical_exp_gene <- clinical_exp_gene %>% 
    dplyr::mutate(exp_group = case_when(exp_value > median_exp ~ "exp_high",
                                        exp_value < median_exp ~ "exp_low"))
  
  ## Do Kaplan-Meier survival analysis
  
  # Fit a survival curve
  surv_object <- survfit(Surv(time_years, vital_status_binary)~exp_group, data = clinical_exp_gene)
  
  # Test for significant difference in survival between the two groups
  #surv_diff <- survdiff(Surv(time_years, vital_status_binary)~exp_group, data = clinical_exp_gene)
  
  # Visualize Kaplan-Meier survival analysis
  km_plot <- ggsurvplot(surv_object, data = clinical_exp_gene, 
                        risk.table = TRUE, pval = TRUE, conf.int = TRUE) +
    labs(x = "Time [years]")
  
  return(km_plot)
  
}

### ---------------------------------------------------------------------------------------------------------



