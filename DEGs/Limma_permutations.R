#!/usr/bin/env Rscript
################################################################################
# Limma_permutations.R
# Permutation analysis with internal loop: runs 100 permutations
# Saves individual permutation results for QQ plots
################################################################################

library(data.table)
library(limma)
library(BiocParallel)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript Limma_permutation.R <timepoint> <n_perms>\n  timepoint: 6 or 24\n  n_perms: number of permutations (default 100)")
}

tp       <- args[1]
n_perms  <- as.integer(args[2])

# Configuration
base_dir   <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates"
output_dir <- file.path(base_dir, "DEG", "permutations")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
today <- format(Sys.time(), "%m%d%Y")

# Parallelization
param <- MulticoreParam(10, progressbar = FALSE)

cat("\n=== PERMUTATION ANALYSIS ===\n")
cat("Timepoint: T", tp, "\n")
cat("Number of permutations:", n_perms, "\n")

################################################################################
# Load expression data (ONCE)
################################################################################

cat("\n=== Loading expression data ===\n")

input_file <- file.path(base_dir, "GxP_Filtered.RData")
load(input_file)
# Provides: filtered_data (list with $dge), anno, c2

dge <- filtered_data$dge

# Subset to this timepoint and relevant treatments
tp_samples <- rownames(dge$samples[dge$samples$Timepoint == tp &
                                   dge$samples$Treatment %in% c("MBP_500nM", "EtOH"), ])
if (length(tp_samples) == 0) stop(paste("No samples found for timepoint", tp))

dge <- dge[, tp_samples]
cat("Samples at timepoint", tp, ":", ncol(dge), "\n")
cat("Genes:", nrow(dge), "\n")
print(table(dge$samples$Treatment))

# Get unique individuals
unique_individuals <- unique(dge$samples$dbGaP_ID)
n_individuals      <- length(unique_individuals)
cat("Total individuals:", n_individuals, "\n")

################################################################################
# Run permutations
################################################################################

cat("\n=== Running", n_perms, "permutations ===\n")

perm_n_degs  <- numeric(n_perms)
perm_top_pval <- numeric(n_perms)

for (perm in 1:n_perms) {

  if (perm %% 10 == 0) cat("Permutation", perm, "/", n_perms, "\n")

  set.seed(perm)

  # Permute treatment within each individual (preserves paired structure)
  dge_perm <- dge
  for (ind in unique_individuals) {
    ind_rows <- which(dge_perm$samples$dbGaP_ID == ind)
    if (length(ind_rows) == 2 && rbinom(1, 1, 0.5) == 1) {
      dge_perm$samples$Treatment[ind_rows] <- rev(dge_perm$samples$Treatment[ind_rows])
    }
  }

  # Set factor levels and make individual ID R-safe
  dge_perm$samples$Treatment     <- factor(dge_perm$samples$Treatment,
                                           levels = c("EtOH", "MBP_500nM"))
  dge_perm$samples$dbGaP_ID_safe <- droplevels(factor(make.names(dge_perm$samples$dbGaP_ID)))

  # Design: treatment + paired individual
  design_perm <- model.matrix(~ Treatment + dbGaP_ID_safe, data = dge_perm$samples)
  cat("Samples:", ncol(dge_perm), "Design rows:", nrow(design_perm), "Design cols:", ncol(design_perm), "\n")

  # Guard: skip if design rows don't match samples (singular permutation)
  if (nrow(design_perm) != ncol(dge_perm)) {
    cat("Skipping permutation", perm, "- design matrix mismatch\n")
    next
  }

  # Run voom on permuted data
  voom_perm <- suppressWarnings(voom(dge_perm, design_perm, plot = FALSE))
  rownames(voom_perm$weights) <- rownames(voom_perm$E)
  colnames(voom_perm$weights) <- colnames(voom_perm$E)

  # Fit model
  fit <- suppressWarnings(limma::lmFit(voom_perm, design_perm, BPPARAM = param))
  fit <- suppressWarnings(limma::eBayes(fit, robust = FALSE))

  # Extract results
  results <- topTable(fit, coef = "TreatmentMBP_500nM", number = Inf, sort.by = "none")

  # Save individual permutation result
  perm_file <- file.path(output_dir,
                         paste0(today, "_MBP_T", tp, "_PERM", perm, ".txt"))
  write.table(results, perm_file, quote = FALSE, sep = "\t",
              row.names = TRUE, col.names = NA)

  # Store summary stats
  perm_n_degs[perm]   <- sum(results$adj.P.Val < 0.1, na.rm = TRUE)
  perm_top_pval[perm] <- min(results$adj.P.Val, na.rm = TRUE)
}

################################################################################
# Save permutation summary
################################################################################

cat("\n=== Saving permutation summary ===\n")

perm_summary <- data.frame(
  permutation = 1:n_perms,
  n_degs      = perm_n_degs,
  top_pval    = perm_top_pval
)

output_file <- file.path(output_dir,
                         paste0(today, "_MBP_T", tp, "_permutation_summary.txt"))
write.table(perm_summary, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Saved:", output_file, "\n")

################################################################################
# Summary statistics
################################################################################

cat("\n=== PERMUTATION SUMMARY ===\n")
cat("Number of DEGs across permutations:\n")
cat("  Mean:",   mean(perm_n_degs),   "\n")
cat("  Median:", median(perm_n_degs), "\n")
cat("  Min:",    min(perm_n_degs),    "\n")
cat("  Max:",    max(perm_n_degs),    "\n")
cat("  SD:",     sd(perm_n_degs),     "\n")

cat("\nTop p-value across permutations:\n")
cat("  Mean:",   mean(perm_top_pval),   "\n")
cat("  Median:", median(perm_top_pval), "\n")
cat("  Min:",    min(perm_top_pval),    "\n")

# Histogram of null distribution
png(file.path(output_dir, paste0(today, "_MBP_T", tp, "_permutation_hist.png")),
    width = 800, height = 600)
hist(perm_n_degs, breaks = 30,
     main  = paste("MBP vs EtOH T", tp, "- Permutation Null Distribution"),
     xlab  = "Number of DEGs (FDR < 0.1)",
     col   = "lightblue",
     border = "black")
dev.off()

cat("\n=== COMPLETE ===\n")
cat("Individual permutation files saved for QQ plot generation\n")
