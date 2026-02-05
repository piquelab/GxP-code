#!/usr/bin/env Rscript
################################################################################
# Prep_genotypes.R
# Process genotype data: PLINK formatting, LD pruning, PCA
# Date: 02/05/2026
################################################################################

# Configuration
BASE <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026"
VCF <- file.path(BASE, "GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz")
OUT_DIR <- file.path(BASE, "tensor")
dir.create(OUT_DIR, showWarnings = FALSE)

OUT_PREFIX <- file.path(OUT_DIR, "GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr")
PRUNE_PREFIX <- file.path(OUT_DIR, "GxP-eQTL_genotypes_pruned")
PCA_PREFIX <- file.path(OUT_DIR, "genotype_pcs")

################################################################################
# Step 1: Convert VCF to PLINK format
################################################################################

cat("Converting VCF to PLINK format...\n")

plink_cmd <- paste(
  "plink2",
  "--vcf", VCF,
  "--max-alleles 2",
  "--maf 0.1",
  "--set-missing-var-ids '@:#:$r:$a'",
  "--chr 1-22",
  "--make-pgen",
  "--out", OUT_PREFIX
)

system(plink_cmd)
cat("PLINK files created with prefix:", OUT_PREFIX, "\n")

################################################################################
# Step 2: LD pruning
################################################################################

cat("\nPerforming LD pruning...\n")

prune_cmd <- paste(
  "plink2",
  "--pfile", OUT_PREFIX,
  "--indep-pairwise 50 5 0.2",
  "--out", PRUNE_PREFIX
)

system(prune_cmd)
cat("LD pruning complete\n")

################################################################################
# Step 3: PCA on pruned variants
################################################################################

cat("\nCalculating PCs...\n")

pca_cmd <- paste(
  "plink2",
  "--pfile", OUT_PREFIX,
  "--extract", paste0(PRUNE_PREFIX, ".prune.in"),
  "--pca 5",
  "--out", PCA_PREFIX
)

system(pca_cmd)
cat("PCA complete\n")

################################################################################
# Step 4: Format PCs for tensorQTL
################################################################################

cat("\nFormatting PCs for tensorQTL...\n")

# Read eigenvec file
eigenvec_file <- paste0(PCA_PREFIX, ".eigenvec")
pcs <- read.table(eigenvec_file, header = TRUE, comment.char = "")

# First column is "#IID" which R reads as "X.IID"
colnames(pcs)[1] <- "IID"

cat("Read", nrow(pcs), "samples with", ncol(pcs) - 1, "PCs\n")
cat("Sample IDs (first 5):", head(pcs$IID, 5), "\n")

# Transpose: PCs as rows, samples as columns
pcs_data <- pcs[, -1]
pcs_t <- t(pcs_data)
colnames(pcs_t) <- pcs$IID
rownames(pcs_t) <- paste0("genoPC", 1:5)

# Create output with id column
pcs_out <- data.frame(id = rownames(pcs_t), pcs_t, check.names = FALSE)

# Save
output_file <- file.path(OUT_DIR, "genotype_pcs_tensorqtl.txt")
write.table(pcs_out, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nFormatted for tensorQTL!\n")
cat("Saved to:", output_file, "\n")
cat("Dimensions:", nrow(pcs_out), "PCs x", ncol(pcs_out) - 1, "samples\n")
cat("Sample IDs (first 5):", head(colnames(pcs_out)[-1], 5), "\n")

cat("\n=== Complete ===\n")
