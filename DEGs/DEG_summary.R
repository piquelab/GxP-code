#!/usr/bin/env Rscript
################################################################################
# DEG_summary.R
# Purpose: Identify DEGs at 10% FDR and create summary tables
################################################################################

library(dplyr)

# Configuration
deg_dir <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/DEG"
fdr_threshold <- 0.10
today <- format(Sys.time(), "%m%d%Y")

contrasts <- c("MBP")
timepoints <- c(6, 24)

################################################################################
# Part 1: DEG outputs
################################################################################

cat("Processing DEG results...\n")

# Collect all DEGs
all_degs <- list()

for (contrast in contrasts) {
  for (tp in timepoints) {
    # Find file
    pattern <- paste0("_", contrast, "_T", tp, ".txt")
    file <- list.files(deg_dir, pattern = pattern, full.names = TRUE)
    
    if (length(file) == 0) {
      cat("Warning: No file found for", contrast, "T", tp, "\n")
      next
    }
    
    # Read results
    res <- read.table(file[1], header = TRUE, sep = "\t", row.names = 1)
    
    # Filter for FDR < 10%
    degs <- res[res$adj.P.Val < fdr_threshold, ]
    
    if (nrow(degs) > 0) {
      degs$Gene <- rownames(degs)
      degs$Contrast <- contrast
      degs$Timepoint <- paste0("T", tp)
      
      # Select columns
      degs <- degs[, c("Gene", "Contrast", "Timepoint", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
      
      all_degs[[paste(contrast, tp, sep = "_")]] <- degs
    }
    
    cat(contrast, "T", tp, ":", nrow(degs), "DEGs\n")
  }
}

# Combine all DEGs
all_degs_df <- do.call(rbind, all_degs)
rownames(all_degs_df) <- NULL

# Save all DEGs file
output_file <- file.path(deg_dir, paste0(today, "_AllDEGs_FDR10.txt"))
write.table(all_degs_df, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Saved all DEGs to:", output_file, "\n")

# Create summary table for DEG
summary_deg <- data.frame(
  BPA_T6 = 0,
  BPA_T24 = 0,
  MBP_T6 = 0,
  MBP_T24 = 0
)

for (contrast in contrasts) {
  for (tp in timepoints) {
    pattern <- paste0("_", contrast, "_T", tp, ".txt")
    file <- list.files(deg_dir, pattern = pattern, full.names = TRUE)
    
    if (length(file) > 0) {
      res <- read.table(file[1], header = TRUE, sep = "\t", row.names = 1)
      n_degs <- sum(res$adj.P.Val < fdr_threshold, na.rm = TRUE)
      
      col_name <- paste0(contrast, "_T", tp)
      summary_deg[1, col_name] <- n_degs
    }
  }
}

# Save summary table
output_file <- file.path(deg_dir, paste0(today, "_DEG_Summary_FDR10.txt"))
write.table(summary_deg, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Saved DEG summary to:", output_file, "\n")

cat("\n=== Complete ===\n")
cat("Total DEGs identified (DEG):", nrow(all_degs_df), "\n")
