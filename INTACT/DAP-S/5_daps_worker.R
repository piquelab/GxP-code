#!/usr/bin/env Rscript
# 5_daps_worker.R
# Run DAP-S fine-mapping for a chunk of genes in one condition.
#
# Usage: Rscript 5_daps_worker.R {condition} {chunk_file}
#
# Inputs:
#   - torus_input/split/{condition}.eQTL.txt.gz  — z-scores
#   - torus_output/{condition}_dump_prior/       — TORUS priors
#   - pgen files                                 — for LD computation
#
# Output:
#   - daps_output/{condition}/{condition}_chunk{N}_pips.txt.gz
#     columns: gene  SC_ID  SNP_ID  PIP

library(data.table)
library(DAPS)

###############################################################################
# CONFIG
###############################################################################

INDIR    <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_input"
TORUSDIR <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_output"
OUTDIR   <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/daps_output"
PGEN     <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/tensor/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"
PLINK2   <- "/wsu/el7/pre-compiled/plink/2.0/plink2"
N        <- 105L   # sample size
TMPDIR   <- tempdir()

###############################################################################
# Parse arguments
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 5_daps_worker.R {condition} {chunk_file}")
CONDITION  <- args[1]
CHUNK_FILE <- args[2]

cat("Condition:", CONDITION, "\n")
cat("Chunk file:", CHUNK_FILE, "\n")

dir.create(file.path(OUTDIR, CONDITION), showWarnings = FALSE, recursive = TRUE)

chunk_id <- sub(".*splitGene", "", sub("\\.txt$", "", basename(CHUNK_FILE)))
outfile  <- file.path(OUTDIR, CONDITION,
                      paste0(CONDITION, "_chunk", chunk_id, "_pips.txt.gz"))

if (file.exists(outfile)) {
  cat("Output already exists, skipping:", outfile, "\n")
  quit(status = 0)
}

###############################################################################
# Load eQTL z-scores for this condition (full file, filter per gene in loop)
###############################################################################

cat("Loading eQTL z-scores...\n")
eqtl <- fread(file.path(INDIR, "split", paste0(CONDITION, ".eQTL.txt.gz")),
              header = TRUE, fill = TRUE)
setnames(eqtl, c("SNP", "gene", "beta", "t.stat", "p.value"))
setkey(eqtl, gene)

genes <- readLines(CHUNK_FILE)
cat("Genes in chunk:", length(genes), "\n\n")

###############################################################################
# Process each gene
###############################################################################

results <- list()

for (gene in genes) {

  cat("  Gene:", gene, "")

  # --- z-scores for this gene ---
  zsub <- eqtl[gene]
  if (nrow(zsub) == 0) { cat("no eQTL data, skipping\n"); next }

  snps <- zsub$SNP
  z    <- zsub$t.stat
  names(z) <- snps

  # --- TORUS prior ---
  prior_file <- file.path(TORUSDIR,
                           paste0(CONDITION, "_dump_prior"),
                           paste0(gene, ".prior"))
  if (!file.exists(prior_file)) { cat("no prior file, skipping\n"); next }

  prior_df <- fread(prior_file, header = FALSE,
                    col.names = c("SNP", "prior"))
  # Align prior to SNP order in z-scores; default to uniform for any missing
  prior_vec <- rep(1 / length(snps), length(snps))
  names(prior_vec) <- snps
  matched <- intersect(snps, prior_df$SNP)
  if (length(matched) > 0) {
    prior_vec[matched] <- prior_df$prior[match(matched, prior_df$SNP)]
  }

  # --- Extract genotypes with plink2 ---
  snp_file <- file.path(TMPDIR, paste0(gene, "_snps.txt"))
  geno_prefix <- file.path(TMPDIR, gene)
  writeLines(snps, snp_file)

  plink_cmd <- paste(
    PLINK2,
    "--pgen", paste0(PGEN, ".pgen"),
    "--pvar", paste0(PGEN, ".pvar"),
    "--psam", paste0(PGEN, ".psam"),
    "--extract", snp_file,
    "--export", "A",          # dosage matrix: one row per sample
    "--out", geno_prefix,
    "--no-psam-pheno",
    "--memory 4000"           # cap RAM per plink2 call
  )
  system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  geno_file <- paste0(geno_prefix, ".raw")
  if (!file.exists(geno_file)) { cat("plink2 failed, skipping\n"); next }

  geno <- fread(geno_file, header = TRUE, fill = TRUE)
  # Remove metadata columns (FID IID PAT MAT SEX PHENOTYPE)
  geno_mat <- as.matrix(geno[, 7:ncol(geno)])

  # plink2 appends _REF to column names — strip to get SNP IDs
  colnames(geno_mat) <- sub("_[^_]+$", "", colnames(geno_mat))

  # Keep only SNPs present in both z-scores and genotype matrix
  common <- intersect(snps, colnames(geno_mat))
  if (length(common) < 2) { cat("too few common SNPs, skipping\n"); next }

  geno_mat  <- geno_mat[, common, drop = FALSE]
  z_common  <- z[common]
  prior_common <- prior_vec[common]

  # Replace any NA dosages with column mean
  geno_mat <- apply(geno_mat, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE); x
  })

  # --- LD matrix ---
  R <- tryCatch(cor(geno_mat), error = function(e) NULL)
  if (is.null(R) || any(is.na(R))) { cat("LD failed, skipping\n"); next }

  # Regularize diagonal to ensure positive definiteness
  diag(R) <- 1.0

  # --- DAP-S ---
  rst <- tryCatch(
    daps_rss(z = z_common, R = R, n = N,
             prior = prior_common,
             estimate_residual_variance = TRUE),
    error = function(e) { cat("DAP-S error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(rst)) next

  # --- Extract PIPs ---
  enloc <- tryCatch(get_enloc(rst), error = function(e) NULL)
  if (is.null(enloc) || nrow(enloc) == 0) {
    # Fall back to raw PIPs if no signal clusters
    enloc <- data.frame(SC_ID = NA, SNP_ID = seq_along(common),
                        PIP = as.numeric(rst$pip))
  }
  # get_enloc returns integer indices for SNP_ID — map back to variant IDs
  enloc$SNP_ID <- common[as.integer(enloc$SNP_ID)]
  enloc$gene <- gene
  results[[gene]] <- enloc[, c("gene", "SC_ID", "SNP_ID", "PIP")]

  # Clean up temp files
  suppressWarnings(file.remove(c(snp_file, geno_file,
                                paste0(geno_prefix, ".log"))))

  cat("-> PIPs for", nrow(enloc), "variants\n")
}

###############################################################################
# Write output
###############################################################################

if (length(results) > 0) {
  out <- rbindlist(results, fill = TRUE)
  fwrite(out, outfile, sep = "\t", compress = "gzip")
  cat("\nWritten:", outfile, "\n")
  cat("Genes processed:", length(results), "/", length(genes), "\n")
} else {
  cat("\nNo results for this chunk.\n")
}
