#!/usr/bin/env Rscript
################################################################################
# Dream.R
# Differential expression analysis using limma on quantile-normalized data
# Uses the same quantile-normalized expression as tensorQTL
################################################################################

library(data.table)
library(limma)
library(variancePartition)
library(BiocParallel)

# Get command line argument for timepoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript Limma.R T6 or Rscript Dream.R T24")
}

tp_arg <- args[1]
tp <- gsub("T", "", tp_arg)

# Configuration
base_dir <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_022526"
output_dir <- file.path(base_dir, "DEG")
dir.create(output_dir, showWarnings = FALSE)
today <- format(Sys.time(), "%m%d%Y")

# Parallelization
param <- MulticoreParam(10, progressbar = FALSE)

cat("\n=== Loading quantile-normalized expression data ===\n")
cat("Timepoint: T", tp, "\n")

################################################################################
# Load expression data from qnorm BED files
################################################################################

conditions <- list(
  BPA_100nM = "BPA",
  H2O = "H2O", 
  MBP_500nM = "MBP",
  EtOH = "EtOH"
)

expr_list <- list()

# Load metadata first to get dbGaP_ID to Sample_ID mapping
input_file <- file.path(base_dir, "GxP_Filtered_02182026.RData")
load(input_file)
meta <- c2[c2$Timepoint == tp, ]

for (cond_name in names(conditions)) {
  cond <- conditions[[cond_name]]
  bed_file <- file.path(base_dir, paste0("GxP-eQTL_", cond, "_T", tp, "_qnorm.bed.gz"))
  
  if (!file.exists(bed_file)) {
    stop(paste("BED file not found:", bed_file))
  }
  
  cat("Loading:", basename(bed_file), "\n")
  bed_data <- fread(cmd = paste("zcat", bed_file), header = TRUE)
  
  # Extract expression matrix (skip first 4 columns: chr, start, end, gene_id)
  gene_ids <- bed_data[[4]]
  expr_matrix <- as.matrix(bed_data[, 5:ncol(bed_data)])
  rownames(expr_matrix) <- gene_ids
  
  # Column names are dbGaP_IDs
  dbgap_ids <- colnames(bed_data)[5:ncol(bed_data)]
  
  # Map dbGaP_IDs to Sample_IDs for this treatment and timepoint
  meta_subset <- meta[meta$Treatment == cond_name, ]
  
  # Create mapping: dbGaP_ID -> Sample_ID
  dbgap_to_sample <- setNames(rownames(meta_subset), meta_subset$dbGaP_ID)
  
  # Rename columns to Sample_IDs
  sample_ids <- dbgap_to_sample[dbgap_ids]
  
  # Check for missing mappings
  if (any(is.na(sample_ids))) {
    warning(paste("Some dbGaP_IDs in", cond_name, "not found in metadata"))
    # Keep original dbGaP_ID if no match
    sample_ids[is.na(sample_ids)] <- dbgap_ids[is.na(sample_ids)]
  }
  
  colnames(expr_matrix) <- sample_ids
  
  expr_list[[cond_name]] <- expr_matrix
  
  cat("  Genes:", nrow(expr_matrix), "\n")
  cat("  Samples:", ncol(expr_matrix), "\n")
}

################################################################################
# Check that all conditions have the same genes
################################################################################

cat("\n=== Validating data ===\n")

# Check genes match
all_genes <- lapply(expr_list, rownames)
if (length(unique(lapply(all_genes, length))) != 1) {
  stop("Different number of genes across conditions")
}
if (!all(sapply(all_genes[-1], function(x) identical(x, all_genes[[1]])))) {
  stop("Gene IDs don't match across conditions")
}

cat("All conditions have", nrow(expr_list[[1]]), "genes\n")

################################################################################
# Combine into single expression matrix
################################################################################

cat("\n=== Combining expression data ===\n")

# Combine all expression data
expr_combined <- do.call(cbind, expr_list)
cat("Combined matrix:", nrow(expr_combined), "genes x", ncol(expr_combined), "samples\n")

################################################################################
# Create sample metadata
################################################################################

cat("\n=== Creating sample metadata ===\n")

# Get Sample_IDs from expression matrix (column names)
sample_ids <- colnames(expr_combined)

# Create metadata from the loaded c2
samples_df <- meta[match(sample_ids, rownames(meta)), ]

# Make sure we have all samples
if (any(is.na(samples_df$Treatment))) {
  stop("Some samples not found in metadata")
}

# Make dbGaP_ID R-safe names for formula
samples_df$dbGaP_ID_safe <- make.names(samples_df$dbGaP_ID)

cat("Sample metadata:\n")
print(table(samples_df$Treatment))


################################################################################
# Run separate analyses for BPA/H2O and MBP/EtOH
################################################################################

# Analysis 1: BPA vs H2O
cat("\n=== ANALYSIS 1: BPA vs H2O ===\n")

samples_bpa           <- samples_df[samples_df$Treatment %in% c("BPA_100nM", "H2O"), ]
expr_bpa              <- expr_combined[, rownames(samples_bpa)]
samples_bpa$Treatment <- factor(samples_bpa$Treatment, levels = c("H2O", "BPA_100nM"))

# EList with uniform weights — appropriate for qnorm data
elist_bpa          <- new("EList")
elist_bpa$E        <- expr_bpa
elist_bpa$weights  <- matrix(1, nrow = nrow(expr_bpa), ncol = ncol(expr_bpa))
elist_bpa$targets  <- samples_bpa

f_bpa <- ~ (1 | Treatment) + (1 | dbGaP_ID)

L_bpa <- makeContrastsDream(f_bpa, samples_bpa, 
  contrasts = c(BPA = "BPA_100nM - H2O"))

fit_bpa <- dream(elist_bpa, f_bpa, samples_bpa, L = L_bpa, ddf = "Satterthwaite", BPPARAM = param)
fit_bpa <- eBayes(fit_bpa)

bpa_results <- topTable(fit_bpa, coef = "BPA", number = Inf, sort.by = "none")

output_file <- file.path(output_dir, paste0(today, "_BPA_T", tp, "_random_dream.txt"))
write.table(bpa_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cat("Saved:", output_file, "\n")
cat("Significant DEGs (FDR < 0.1):", sum(bpa_results$adj.P.Val < 0.1), "\n")

################################################################################
# Analysis 2: MBP vs EtOH
################################################################################

cat("\n=== ANALYSIS 2: MBP vs EtOH ===\n")

samples_mbp           <- samples_df[samples_df$Treatment %in% c("MBP_500nM", "EtOH"), ]
expr_mbp              <- expr_combined[, rownames(samples_mbp)]
samples_mbp$Treatment <- factor(samples_mbp$Treatment, levels = c("EtOH", "MBP_500nM"))

elist_mbp         <- new("EList")
elist_mbp$E       <- expr_mbp
elist_mbp$weights <- matrix(1, nrow = nrow(expr_mbp), ncol = ncol(expr_mbp))
elist_mbp$targets <- samples_mbp

f_mbp <- ~ (1 | Treatment) + (1 | dbGaP_ID)

L_mbp <- makeContrastsDream(f_mbp, samples_mbp, 
  contrasts = c(MBP = "MBP_500nM - EtOH"))

fit_mbp <- dream(elist_mbp, f_mbp, samples_mbp, L = L_mbp, ddf = "Satterthwaite", BPPARAM = param)
fit_mbp <- eBayes(fit_mbp)

mbp_results <- topTable(fit_mbp, coef = "MBP", number = Inf, sort.by = "none")

output_file <- file.path(output_dir, paste0(today, "_MBP_T", tp, "_random_dream.txt"))
write.table(mbp_results, output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cat("Saved:", output_file, "\n")
cat("Significant DEGs (FDR < 0.1):", sum(mbp_results$adj.P.Val < 0.1), "\n")

cat("\n=== COMPLETE ===\n")
cat("\n=== COMPLETE ===\n")
