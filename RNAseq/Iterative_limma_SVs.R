#!/usr/bin/env Rscript
################################################################################
# Iterative_limma_SVs.R
# Purpose: Iteratively include SVs in voom normalization and limma analysis
# Date: 02/05/2026
################################################################################

library(limma)
library(edgeR)
library(BiocParallel)

# Configuration
input_file <- "GxP_Filtered_02052026.RData"
output_dir <- "DEG_SVs"
dir.create(output_dir, showWarnings = FALSE)
today <- format(Sys.time(), "%m%d%Y")
timepoints <- c(6, 24)
max_svs <- 25

# Parallelization
param <- MulticoreParam(10, progressbar = FALSE)

# Load filtered data
cat("Loading filtered data...\n")
load(input_file)

# Process each timepoint
for (tp in timepoints) {
  cat("\n=== Timepoint", tp, "===\n")
  
  # Load DGE object (pre-voom)
  tp_data <- filtered_data[[paste0("tp", tp)]]
  dge <- tp_data$dge
  samples <- dge$samples
  
  # Load SVs
  sv_file <- file.path("Filtered_SVs", paste0("SV_tp", tp, ".RData"))
  cat("Loading SVs from:", sv_file, "\n")
  load(sv_file)
  
  # Add SVs to sample metadata
  samples <- cbind(samples, sv_matrix)
  
  # Set factor levels
  samples$Treatment <- factor(samples$Treatment, levels = c("H2O", "EtOH", "BPA_100nM", "MBP_500nM"))
  samples$dbGaP_ID <- make.names(samples$dbGaP_ID)
  
  # Iterate through SV combinations (0 to 25)
  for (n_sv in 0:max_svs) {
    cat("\n--- Including", n_sv, "SVs ---\n")
    
    # Build formula
    if (n_sv == 0) {
      f <- as.formula("~ Treatment + dbGaP_ID + trimmed_dClean.dFastq")
      sv_label <- "SV0"
    } else {
      sv_terms <- paste0("SV", 1:n_sv, collapse = " + ")
      f <- as.formula(paste("~ Treatment + dbGaP_ID + trimmed_dClean.dFastq +", sv_terms))
      sv_label <- paste0("SV1-", n_sv)
    }
    
    cat("Formula:", deparse(f), "\n")
    
    # Create design matrix
    design <- model.matrix(f, data = samples)
    
    # Voom transformation with design
    voom_data <- voom(dge, design = design, plot = FALSE)
    
    # Fit limma model
    fit <- limma::lmFit(voom_data, design, BPPARAM = param)
    fit <- limma::eBayes(fit, robust = FALSE)
    
    # Extract BPA vs H2O
    bpa_results <- topTable(fit, coef = "TreatmentBPA_100nM", number = Inf, sort.by = "none")
    output_file <- file.path(output_dir, paste0(today, "_BPA_T", tp, "_", sv_label, ".txt"))
    write.table(bpa_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
    
    # Extract EtOH vs H2O
    etoh_results <- topTable(fit, coef = "TreatmentEtOH", number = Inf, sort.by = "none")
    output_file <- file.path(output_dir, paste0(today, "_CRL_T", tp, "_", sv_label, ".txt"))
    write.table(etoh_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
    
    # Extract MBP vs EtOH (contrast)
    colnames(fit$coefficients)[1] <- "Intercept"
    colnames(fit$design)[1] <- "Intercept"
    contrast_matrix <- makeContrasts(
      MBP = TreatmentMBP_500nM - TreatmentEtOH,
      levels = colnames(fit$coefficients)
    )
    
    fit_contrast <- contrasts.fit(fit, contrast_matrix)
    fit_contrast <- eBayes(fit_contrast, robust = FALSE)
    
    mbp_results <- topTable(fit_contrast, coef = "MBP", number = Inf, sort.by = "none")
    output_file <- file.path(output_dir, paste0(today, "_MBP_T", tp, "_", sv_label, ".txt"))
    write.table(mbp_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
    
    cat("Saved results for", sv_label, "\n")
  }
}

cat("\n=== Complete ===\n")
cat("Total files created:", 2 * (max_svs + 1) * 3, "\n")
