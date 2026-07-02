#!/usr/bin/env Rscript
################################################################################
# Limma.R
# Differential expression analysis using limma on voom-normalized data
# Uses GxP_Filtered.RData output from GxP_eQTL_prep.R
# Voom is run per-timepoint contrast for correct mean-variance estimation
################################################################################

library(data.table)
library(limma)
library(BiocParallel)

# Get command line argument for timepoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript Limma.R T6 or Rscript Limma.R T24")
}

tp_arg <- args[1]
tp <- gsub("T", "", tp_arg)

# Configuration
base_dir   <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates"
output_dir <- file.path(base_dir, "DEG")
dir.create(output_dir, showWarnings = FALSE)
today <- format(Sys.time(), "%m%d%Y")

# Parallelization
param <- MulticoreParam(10, progressbar = FALSE)

################################################################################
# Load filtered data
################################################################################

cat("\n=== Loading filtered expression data ===\n")
cat("Timepoint: T", tp, "\n")

input_file <- file.path(base_dir, "GxP_Filtered.RData")
load(input_file)
# Provides: filtered_data (list with $dge), anno, c2

dge <- filtered_data$dge

# Subset to this timepoint and relevant treatments
tp_samples  <- rownames(dge$samples[dge$samples$Timepoint == tp &
                                    dge$samples$Treatment %in% c("MBP_500nM", "EtOH"), ])
if (length(tp_samples) == 0) stop(paste("No samples found for timepoint", tp))

dge <- dge[, tp_samples]
cat("Samples at timepoint", tp, ":", ncol(dge), "\n")
cat("Genes:", nrow(dge), "\n")
print(table(dge$samples$Treatment))

################################################################################
# Analysis: MBP vs EtOH
################################################################################

cat("\n=== ANALYSIS: MBP vs EtOH ===\n")

# Set reference level and make individual ID R-safe
dge$samples$Treatment     <- factor(dge$samples$Treatment, levels = c("EtOH", "MBP_500nM"))
dge$samples$dbGaP_ID_safe <- make.names(dge$samples$dbGaP_ID)

# Design: treatment + paired individual + batch effect
design_mbp <- model.matrix(~ Treatment + dbGaP_ID_safe + trimmed_dClean.dFastq, data = dge$samples)
cat("Design matrix:", nrow(design_mbp), "samples x", ncol(design_mbp), "coefficients\n")

# Run voom on this contrast's samples only
cat("Running voom...\n")
voom_data <- voom(dge, design_mbp, plot = FALSE)

# Ensure dimnames are set on weights
rownames(voom_data$weights) <- rownames(voom_data$E)
colnames(voom_data$weights) <- colnames(voom_data$E)

# Fit model
fit_mbp <- limma::lmFit(voom_data, design_mbp, BPPARAM = param)
fit_mbp <- limma::eBayes(fit_mbp, robust = FALSE)

# Extract results
mbp_results <- topTable(fit_mbp, coef = "TreatmentMBP_500nM", number = Inf, sort.by = "none")
output_file <- file.path(output_dir, paste0(today, "_MBP_T", tp, ".txt"))
write.table(mbp_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cat("Saved:", output_file, "\n")
cat("Significant DEGs (FDR < 0.05):", sum(mbp_results$adj.P.Val < 0.05), "\n")
cat("Significant DEGs (FDR < 0.10):", sum(mbp_results$adj.P.Val < 0.10), "\n")

cat("\n=== COMPLETE ===\n")
