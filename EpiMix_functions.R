## All below functions are from EpiMix:
## https://github.com/gevaertlab/EpiMix
## Zheng Y, Jun J, Gevaert O (2023). EpiMix: 
## EpiMix: an integrative tool for the population-level 
## analysis of DNA methylation. R package version 1.2.4.
## with associated publication:
## Zheng et al. 2023 EpiMix is an integrative tool for epigenomic subtyping 
## using DNA methylation, Cell Reports Methods, doi: 10.1016/j.crmeth.2023.100515

#' The Preprocess_DNAMethylation function
#'
#' @description Preprocess DNA methylation data from the GEO database.
#' @details
#' The data preprocessing pipeline includes:
#' (1) eliminating samples and genes with too many NAs, imputing NAs.
#' (2) (optional) mapping the column names of the DNA methylation data to the actual sample names based on the information from 'sample.map'.
#' (3) (optional) removing CpG probes on the sex chromosomes or the user-defined chromosomes.
#' (4) (optional) doing Batch correction.
#' If both sample.info and group.1 and group.2 information are provided, the function will perform missing value estimation and batch correction on group.1 and group.2 separately. This will ensure that the true difference between group.1 and group.2 will not be obscured by missing value estimation and batch correction.
#' @param methylation.data matrix of DNA methylation data with CpG in rows and sample names in columns.
#' @param met.platform character string indicating the type of the Illumina Infinium BeadChip for collecting the methylation data. Should be either 'HM450' or 'EPIC'. Default: 'EPIC'
#' @param genome character string indicating the genome build version for retrieving the probe annotation. Should be either 'hg19' or 'hg38'. Default: 'hg38'.
#' @param sample.info dataframe that maps each sample to a study group. Should contain two columns: the first column (named: 'primary') indicating the sample names, and the second column (named: 'sample.type') indicating which study group each sample belongs to (e.g., “Experiment” vs. “Control”,  “Cancer” vs. “Normal”). Sample names in the 'primary' column must coincide with the column names of the methylation.data. Please see details for more information. Default: NULL.
#' @param group.1 character vector indicating the name(s) for the experiment group. The values must coincide with the values in the 'sample.type' of the sample.info dataframe.Please see details for more information. Default: NULL.
#' @param group.2 character vector indicating the names(s) for the control group. The values must coincide with the values in the 'sample.type' of the sample.info dataframe. Please see details for more information. Default: NULL.
#' @param sample.map dataframe for mapping the GEO accession ID (column names) to the actual sample names. Can be the output from the GEO_getSampleMap function. Default: NULL.
#' @param rm.chr character vector indicating the probes on which chromosomes to be removed. Default: 'chrX', 'chrY'.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default: 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default: 0.1.
#' @param doBatchCorrection logical indicating whether to perform batch correction. If TRUE, the batch data need to be provided.
#' @param BatchData dataframe with batch information. Should contain two columns: the first column indicating the actual sample names, the second column indicating the batch. Users are expected to retrieve the batch information from the GEO on their own, but this can also be done using the GEO_getSampleInfo function with the 'group.column' as the column indicating the batch for each sample. Defualt': NULL.
#' @param batch.correction.method character string indicating the method that will be used for batch correction. Should be either 'Seurat' or 'Combat'. Default: 'Seurat'.
#' @param cores number of CPU cores to be used for batch effect correction. Defaut: 1.
#' @return DNA methylation data matrix with probes in rows and samples in columns.
#' @import doSNOW
#' @import doParallel
#' @export
#' @keywords preprocess
#' @examples
#' {
#' data(MET.data)
#' data(LUAD.sample.annotation)
#'
#' Preprocessed_Data <- Preprocess_DNAMethylation(MET.data,
#'                                                    met.platform = 'HM450',
#'                                                    sample.info = LUAD.sample.annotation,
#'                                                    group.1 = 'Cancer',
#'                                                    group.2 = 'Normal')
#'
#' }

Preprocess_DNAMethylation <- function(methylation.data, met.platform = "EPIC",
                                      genome = "hg38", sample.info = NULL, group.1 = NULL, group.2 = NULL, sample.map = NULL,
                                      rm.chr = c("chrX", "chrY"), MissingValueThresholdGene = 0.2, MissingValueThresholdSample = 0.2,
                                      doBatchCorrection = FALSE, BatchData = NULL, batch.correction.method = "Seurat",
                                      cores = 1) {
  
  # check data
  if (nrow(methylation.data) <= 1 | ncol(methylation.data) <= 1) {
    stop("methylation.data is empty!\n Please check whether the actual data are saved as supplementary files in GEO.\n")
  }
  
  ### Step 1: convert column names to the actual patient names
  if (!is.null(sample.map)) {
    cat("Mapping column names of the DNA methylation data to the actual sample names...\n")
    sample.map <- sample.map[sample.map$assay == "DNA methylation", ]
    overlapSamples <- intersect(sample.map$colnames, colnames(methylation.data))
    if (length(overlapSamples) > 0) {
      methylation.data <- methylation.data[, overlapSamples, drop = FALSE]
      presentSamples <- match(colnames(methylation.data), sample.map$colnames)
      sampleNames <- sample.map$primary
      sampleNames <- sampleNames[presentSamples]
      colnames(methylation.data) <- sampleNames
    } else {
      warning("No overlap samples were found between the sample map and the column names of the methylation data! No mapping of column names was performed!")
    }
  }
  
  ### Step 2: Filter out SNP probes
  cat("\tFiltering out SNP probes...\n")
  GoodProbes <- rownames(methylation.data)[startsWith(rownames(methylation.data),
                                                      "cg")]
  NrProbesToRemove <- length(rownames(methylation.data)) - length(GoodProbes)
  methylation.data <- methylation.data[GoodProbes, ]
  cat("\tRemoved", NrProbesToRemove, "probes with SNPs.\n")
  
  ### Step 3:Filter out CpG probes in the user-specified chromosomes
  cat("\tFetching probe annotation for", met.platform, "\n")
  ProbeAnnotation <- EpiMix_getInfiniumAnnotation(plat = met.platform, genome = genome)
  ProbeAnnotation <- convertAnnotToDF(ProbeAnnotation)
  ProbesToRemove <- ProbeAnnotation$probeID[ProbeAnnotation$CpG_chrm %in% rm.chr]
  methylation.data <- methylation.data[-which(rownames(methylation.data) %in% ProbesToRemove),
                                       , drop = FALSE]
  cat("Removing", length(which(rownames(methylation.data) %in% ProbesToRemove)),
      "CpG probes on", paste0(rm.chr, ","), "\n")
  
  ### Step 4: Split the experiment and the control group
  MET_Experiment <- MET_Control <- NULL
  if (!is.null(sample.info) & !is.null(group.1) & !is.null(group.2)) {
    overlapSamples <- intersect(sample.info$primary, colnames(methylation.data))
    methylation.data <- methylation.data[, overlapSamples, drop = FALSE]
    sample.info <- sample.info[sample.info$primary %in% overlapSamples, ]
    cat("Found", length(overlapSamples), "samples with sample information.\n")
    Samples_Experiment <- sample.info[sample.info$sample.type %in% group.1, "primary"]
    Samples_Control <- sample.info[sample.info$sample.type %in% group.2, "primary"]
    MET_Experiment <- methylation.data[, Samples_Experiment, drop = FALSE]
    MET_Control <- methylation.data[, Samples_Control, drop = FALSE]
    if (ncol(MET_Experiment) == 0)
      stop("Cannot find methylation data with samples in the group.1 ! The sample names must overlap with the column names of the methylation data.")
    if (ncol(MET_Control) == 0)
      stop("Cannot find methylation data with samples in the group.2 ! The sample names must overlap with the column names of the methylation data.")
    cat("There are", ncol(MET_Experiment), "samples in the", paste0(group.1,
                                                                    collapse = " and "), "group, and", ncol(MET_Control), "samples in the",
        paste0(group.2, collapse = " and "), "group.\n")
  } else {
    warning("sample.info or group.1 and group.2 information is not provided, the funciton will perform the missing value estimation and the batch correction on the entire dataset.")
    MET_Experiment <- methylation.data
  }
  rm(methylation.data)
  gc()
  
  ### Step 4: Missing value estimation
  if (sum(is.na(MET_Experiment)) > 0) {
    cat("\tMissing value estimation on group.1...\n")
    MET_Experiment <- GEO_EstimateMissingValues_Methylation(MET_Experiment, MissingValueThresholdGene,
                                                            MissingValueThresholdSample)
  }
  if (!is.null(MET_Control) & sum(is.na(MET_Control)) >
      0) {
    cat("\tMissing value estimation on group.2...\n")
    MET_Control <- GEO_EstimateMissingValues_Methylation(MET_Control, MissingValueThresholdGene,
                                                         MissingValueThresholdSample)
  }
  
  ### Step 5: Batch correction
  if (doBatchCorrection) {
    MinInBatch <- 0
    if (batch.correction.method == "Combat") {
      MinInBatch <- 5
      if (cores > 1) {
        # unregister()
        cat("Registering sockets on multiple CPU cores...\n")
        cl <- parallel::makeCluster(cores)
      }
    }
    cat("Performing batch correction on group.1...\n")
    MET_Experiment <- CorrectBatchEffect(MET_Experiment, BatchData, batch.correction.method,
                                         MinInBatch = MinInBatch, featurePerSet = 50000)
    if (!is.null(MET_Control) && ncol(MET_Control) > 0) {
      cat("Performing batch correction on group.2...\n")
      MET_Control <- CorrectBatchEffect(MET_Control, BatchData, batch.correction.method,
                                        MinInBatch = MinInBatch, featurePerSet = 50000)
    }
    
    # Set values <0 to 0 and >1 to 1, because of batch correction
    MET_Experiment[MET_Experiment < 0] <- 0
    MET_Experiment[MET_Experiment > 1] <- 1
    if (!is.null(MET_Control) & ncol(MET_Control) > 0) {
      MET_Control[MET_Control < 0] <- 0
      MET_Control[MET_Control > 1] <- 1
    }
    if (cores > 1 & batch.correction.method == "Combat")
      parallel::stopCluster(cl)
  }
  
  
  ### Step 6: combine MET_Experiment and MET_Control into one matrix
  if (!is.null(MET_Control) && ncol(MET_Control) > 0) {
    overlapProbes <- intersect(rownames(MET_Experiment), rownames(MET_Control))
    cat("Found", length(overlapProbes), "overlapping probes between group.1 and group.2 after preprocessing...\n")
    MET_Experiment <- MET_Experiment[overlapProbes, , drop = FALSE]
    MET_Control <- MET_Control[overlapProbes, , drop = FALSE]
    MET_Data <- cbind(MET_Experiment, MET_Control)
    return(MET_Data)
  } else {
    return(MET_Experiment)
  }
}

#' EpiMix_getInfiniumAnnotation
#' This function gets Infinium probe annotation from the
#' sesameData library. This function is from the EpiMix package
#' https://bioconductor.org/packages/release/bioc/html/EpiMix.html.
#' Zheng Y, Jun J, Gevaert O (2023). EpiMix: 
#' EpiMix: an integrative tool for the population-level analysis of 
#' DNA methylation. R package version 1.1.2.
#' @param plat A character string representing the methylation platform 
#' which can either be HM27, HM450 or EPIC
#' @param genome A character string representing the genome build version
#' which can either be hg19 or hg38
#' @import ExperimentHub
#' @return Probe annotations
#' @keywords internal
EpiMix_getInfiniumAnnotation <- function(plat = "EPIC", genome = "hg38") {
  hubID <- NULL
  if (tolower(genome) == "hg19" & toupper(plat) == "HM27")
    hubID <- "EH3672"
  if (tolower(genome) == "hg38" & toupper(plat) == "HM27")
    hubID <- "EH3673"
  if (tolower(genome) == "hg19" & toupper(plat) == "HM450")
    hubID <- "EH3674"
  if (tolower(genome) == "hg38" & toupper(plat) == "HM450")
    hubID <- "EH3675"
  if (tolower(genome) == "hg19" & toupper(plat) == "EPIC")
    hubID <- "EH3670"
  if (tolower(genome) == "hg38" & toupper(plat) == "EPIC")
    hubID <- "EH3671"
  ProbeAnnotation <- ExperimentHub::ExperimentHub()[[hubID]]
  return(ProbeAnnotation)
}

#' The convertAnnotToDF function
#' @description convert the probe annotation from the GRange object to a dataframe
#' @param annot a GRange object of probe annotation, can be the object returned from the getInfiniumAnnotation function.
#' @return a dataframe with chromosome, beginning and end position, mapped gene information for each CpG probe
#' @keywords internal
#'
convertAnnotToDF <- function(annot) {
  df.annot <- data.frame(CpG_chrm = GenomicRanges::seqnames(annot), CpG_beg = GenomicRanges::start(ranges(annot)),
                         CpG_end = GenomicRanges::end(ranges(annot)), probeID = names(annot), gene = GenomicRanges::mcols(annot)$gene)
  return(df.annot)
}


#' The GEO_EstimateMissingValues_Methylation function
#'
#' Internal. Removes samples and probes with more missing values than the MissingValueThreshold, and imputes remaining missing values using Tibshirani's KNN method.
#' @param MET_Data methylation data or gene expression data matrix.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return the dataset with imputed values and possibly some genes or samples deleted.
#' @keywords internal
#'
GEO_EstimateMissingValues_Methylation <- function(MET_Data, MissingValueThresholdGene = 0.3,
                                                  MissingValueThresholdSample = 0.3) {
  
  # removing probes with too many missing values
  NrMissingsPerGene <- apply(MET_Data, 1, function(x) sum(is.na(x)))/ncol(MET_Data)
  cat("Removing", sum(NrMissingsPerGene > MissingValueThresholdGene), "probes with more than",
      MissingValueThresholdGene * 100, "% missing values.\n")
  if (sum(NrMissingsPerGene > MissingValueThresholdGene) > 0)
    MET_Data <- MET_Data[NrMissingsPerGene < MissingValueThresholdGene, ]
  
  # removing samples with too many missings values
  NrMissingsPerSample <- apply(MET_Data, 2, function(x) sum(is.na(x)))/nrow(MET_Data)
  cat("Removing", sum(NrMissingsPerSample > MissingValueThresholdSample), "samples with more than",
      MissingValueThresholdSample * 100, "% missing values.\n")
  if (sum(NrMissingsPerSample > MissingValueThresholdSample) > 0)
    MET_Data <- MET_Data[, NrMissingsPerSample < MissingValueThresholdSample]
  
  # knn impute using Tibshirani's method
  if (length(colnames(MET_Data)) > 1) {
    k <- 15
    KNNresults <- impute::impute.knn(as.matrix(MET_Data), k)
    MET_Data_KNN <- KNNresults$data
    # cleaning up sample names
    return(MET_Data_KNN)
    
  } else {
    # when only 1 sample,need to make a matrix again
    # MET_Data=as.matrix(MET_Data)
    return(MET_Data)
  }
}
