#!/usr/bin/env Rscript
################################################################################
# Filtered_SVA.R
# Purpose: Run surrogate variable analysis on filtered voom-normalized data
# Date: 02/05/2026
################################################################################

library(sva)

# Configuration
input_file <- "GxP_Filtered_02052026.RData"
output_dir <- "Filtered_SVs"
dir.create(output_dir, showWarnings = FALSE)
timepoints <- c(6, 24)

# Load data
cat("Loading data...\n")
load(input_file)

# Process each timepoint
for (tp in timepoints) {
  cat("\n=== Timepoint", tp, "===\n")
  
  tp_data <- filtered_data[[paste0("tp", tp)]]
  voom_data <- tp_data$voom
  samples <- tp_data$dge$samples
  
  # Expression matrix (log2-CPM from voom)
  expr <- voom_data$E
  
  cat("Genes:", nrow(expr), "\n")
  cat("Samples:", ncol(expr), "\n")
  
  # Full model
  mod <- model.matrix(~Treatment + dbGaP_ID + trimmed_dClean.dFastq, data = samples)
  
  # Null model
  mod0 <- model.matrix(~dbGaP_ID + trimmed_dClean.dFastq, data = samples)
  
  # Estimate number of SVs
  cat("Estimating number of surrogate variables...\n")
  n.sv <- num.sv(expr, mod, method = "be")
  cat("Estimated SVs:", n.sv, "\n")
  
  # Run SVA
  cat("Running SVA...\n")
  svobj <- sva(expr, mod, mod0, n.sv = n.sv)
  
  cat("Identified", svobj$n.sv, "surrogate variables\n")
  
  # SV matrix
  sv_matrix <- svobj$sv
  rownames(sv_matrix) <- colnames(expr)
  colnames(sv_matrix) <- paste0("SV", 1:svobj$n.sv)
  
  # Save results
  output_file <- file.path(output_dir, paste0("SV_tp", tp, ".RData"))
  save(svobj, sv_matrix, file = output_file)
  cat("Saved:", output_file, "\n")
  
  output_txt <- file.path(output_dir, paste0("SV_tp", tp, ".txt"))
  write.table(sv_matrix, output_txt, sep = "\t", quote = FALSE, 
              row.names = TRUE, col.names = NA)
  cat("Saved:", output_txt, "\n")
}

cat("\n=== Complete ===\n")
