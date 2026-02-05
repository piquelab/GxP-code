#!/usr/bin/env Rscript
################################################################################
# Limma.R
# Purpose: Differential expression analysis using limma
# Date: 02/05/2026
################################################################################

library(limma)
library(edgeR)
library(BiocParallel)

# Get command line argument for timepoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript Limma.R T6 or Rscript Limma.R T24")
}

tp_arg <- args[1]
tp <- as.numeric(gsub("T", "", tp_arg))

# Configuration
input_file <- "GxP_Filtered_02052026.RData"
output_dir <- "DEG"
dir.create(output_dir, showWarnings = FALSE)
today <- format(Sys.time(), "%m%d%Y")

# Parallelization
param <- MulticoreParam(10, progressbar = FALSE)

# Load data
cat("Loading data...\n")
load(input_file)

cat("\n=== Timepoint", tp, "===\n")

# Get data for this timepoint
tp_data <- filtered_data[[paste0("tp", tp)]]
voom_data <- tp_data$voom
samples <- tp_data$dge$samples

cat("Genes:", nrow(voom_data$E), "\n")
cat("Samples:", ncol(voom_data$E), "\n")

# Set treatment levels (H2O as reference/intercept)
samples$Treatment <- factor(samples$Treatment, levels = c("H2O", "EtOH", "BPA_100nM", "MBP_500nM"))
samples$dbGaP_ID <- make.names(samples$dbGaP_ID)

# Create design matrix with intercept
f <- as.formula("~ Treatment + dbGaP_ID + trimmed_dClean.dFastq")
design <- model.matrix(f, data = samples)

# Fit model
cat("Fitting limma model...\n")
fit <- limma::lmFit(voom_data, design, BPPARAM = param)
fit <- limma::eBayes(fit, robust = FALSE)

# Save BPA vs H2O
cat("Extracting BPA vs H2O...\n")
coef <- "TreatmentBPA_100nM"
bpa_results <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
output_file <- file.path(output_dir, paste0(today, "_BPA_T", tp, ".txt"))
write.table(bpa_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cat("Saved:", output_file, "\n")

# Save EtOH vs H2O (control comparison)
cat("Extracting EtOH vs H2O...\n")
coef <- "TreatmentEtOH"
etoh_results <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
output_file <- file.path(output_dir, paste0(today, "_CRL_T", tp, ".txt"))
write.table(etoh_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cat("Saved:", output_file, "\n")

# MBP vs EtOH contrast
cat("Extracting MBP vs EtOH...\n")

# Rename the intercept
colnames(fit$coefficients)[1] <- "Intercept"
colnames(fit$design)[1] <- "Intercept"

# Use colnames from the fitted model
contrast_matrix <- makeContrasts(
  MBP = TreatmentMBP_500nM - TreatmentEtOH,
  levels = colnames(fit$coefficients)
)

fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_contrast <- eBayes(fit_contrast, robust = FALSE)

mbp_results <- topTable(fit_contrast, coef = "MBP", number = Inf, sort.by = "none")
output_file <- file.path(output_dir, paste0(today, "_MBP_T", tp, ".txt"))
write.table(mbp_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cat("Saved:", output_file, "\n")

cat("\n=== Complete ===\n")
