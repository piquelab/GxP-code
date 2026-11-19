#!/usr/bin/env Rscript
# 5_daps_worker.R
# Run DAP-S fine-mapping for a chunk of genes in one condition.
#
# Usage: Rscript 5_daps_worker.R {condition} {chunk_file}
#
# Inputs:
#   torus_input/split/{condition}.eQTL.txt.gz  — z-scores
#   torus_output/{condition}_dump_prior/        — TORUS priors
#   pgen files                                  — for LD computation
#
# Output:
#   daps_output/{condition}/{condition}_chunk{N}_pips.txt.gz
#   columns: gene  SC_ID  SNP_ID  PIP

library(data.table)
library(DAPS)

# Fix daps_rss bug: prior argument silently dropped
original_body <- deparse(body(daps_rss))
original_body <- sub("prior = NULL, null_wt = NULL,",
                     "prior = prior, null_wt = null_wt,",
                     original_body, fixed = TRUE)
body(daps_rss) <- parse(text = original_body)
environment(daps_rss) <- asNamespace("DAPS")
assignInNamespace("daps_rss", daps_rss, ns = "DAPS")

###############################################################################
# CONFIG
###############################################################################

INDIR    <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/torus_input"
TORUSDIR <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/torus_output/open_chromatin"
OUTDIR   <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/daps_output/open_chromatin"
PGEN     <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/starter_files/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"
PLINK2   <- "/wsu/el7/pre-compiled/plink/2.0/plink2"
N        <- 105L #samples

###############################################################################
# Parse arguments
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 5_daps_worker.R {condition} {chunk_file}")
CONDITION  <- args[1]
CHUNK_FILE <- args[2]

cat("Condition:", CONDITION, "\n")
cat("Chunk file:", CHUNK_FILE, "\n")

# Job-specific temp directory to avoid collisions between parallel jobs
TMPDIR <- file.path(tempdir(),
                    paste0(CONDITION, "_", sub("\\.txt$", "", basename(CHUNK_FILE))))
dir.create(TMPDIR, showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(OUTDIR, CONDITION), showWarnings = FALSE, recursive = TRUE)

chunk_id <- sub(".*splitGene", "", sub("\\.txt$", "", basename(CHUNK_FILE)))
outfile  <- file.path(OUTDIR, CONDITION,
                      paste0(CONDITION, "_chunk", chunk_id, "_pips.txt.gz"))

if (file.exists(outfile)) {
  cat("Output already exists, skipping:", outfile, "\n")
  quit(status = 0)
}

###############################################################################
# Load eQTL z-scores — use setkey for fast per-gene subsetting
###############################################################################

cat("Loading eQTL z-scores...\n")
eqtl <- fread(file.path(INDIR, "split", paste0(CONDITION, ".eQTL.txt.gz")),
              header = TRUE, fill = TRUE)
setnames(eqtl, c("SNP", "gene", "beta", "t.stat", "p.value"))
setkey(eqtl, gene)   # enables fast eqtl[.(gene_id)] lookup

genes <- readLines(CHUNK_FILE)
cat("Genes in chunk:", length(genes), "\n\n")

###############################################################################
# Process each gene
###############################################################################

results <- list()

for (target_gene in genes) {

  cat("  Gene:", target_gene, "")

  # --- z-scores: use .() to avoid column-name scoping ambiguity ---
  zsub <- eqtl[.(target_gene)]
  if (nrow(zsub) == 0) { cat("no eQTL data, skipping\n"); next }

  snps <- zsub$SNP
  z    <- zsub$t.stat
  names(z) <- snps

  # --- TORUS prior ---
  prior_file <- file.path(TORUSDIR,
                           paste0(CONDITION, "_dump_prior"),
                           paste0(target_gene, ".prior"))
  if (!file.exists(prior_file)) { cat("no prior file, skipping\n"); next }

  prior_df <- fread(prior_file, header = FALSE, col.names = c("SNP", "prior"))

  prior_vec <- rep(1 / length(snps), length(snps))
  names(prior_vec) <- snps
  matched <- intersect(snps, prior_df$SNP)
  if (length(matched) > 0) {
    prior_vec[matched] <- prior_df$prior[match(matched, prior_df$SNP)]
  }

  # --- Extract genotypes with plink2 ---
  snp_file    <- file.path(TMPDIR, paste0(target_gene, "_snps.txt"))
  geno_prefix <- file.path(TMPDIR, target_gene)
  writeLines(snps, snp_file)

  plink_cmd <- paste(
    PLINK2,
    "--pgen", paste0(PGEN, ".pgen"),
    "--pvar", paste0(PGEN, ".pvar"),
    "--psam", paste0(PGEN, ".psam"),
    "--extract", snp_file,
    "--export A",
    "--out", geno_prefix,
    "--no-psam-pheno",
    "--memory 12000"
  )
  system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  geno_file <- paste0(geno_prefix, ".raw")
  if (!file.exists(geno_file)) { cat("plink2 failed, skipping\n"); next }

  geno <- fread(geno_file, header = TRUE)
  geno_mat <- as.matrix(geno[, 7:ncol(geno)])

  # Strip plink2 allele suffix from column names (e.g. "1:924320:T:C_C" -> "1:924320:T:C")
  colnames(geno_mat) <- sub("_[^_]+$", "", colnames(geno_mat))

  # Keep only SNPs present in both z-scores and genotype matrix
  common <- intersect(snps, colnames(geno_mat))
  if (length(common) < 2) { cat("too few common SNPs, skipping\n"); next }

  geno_mat     <- geno_mat[, common, drop = FALSE]
  z_common     <- z[common]
  prior_common <- prior_vec[common]

  # Impute missing dosages with column mean
  geno_mat <- apply(geno_mat, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE); x
  })

  # --- LD matrix ---
  R <- tryCatch(cor(geno_mat), error = function(e) NULL)
  if (is.null(R) || any(is.na(R))) { cat("LD failed, skipping\n"); next }
  diag(R) <- 1.0

  # --- DAP-S ---
  rst <- tryCatch(
    daps_rss(z = z_common, R = R, n = N,
             prior = prior_common,
             estimate_residual_variance = TRUE),
    error = function(e) { cat("DAP-S error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(rst)) next

# --- Extract PIPs for ALL variants ---
# Start with full PIP vector for every tested variant
all_pips <- data.frame(
  gene   = target_gene,
  SC_ID  = NA_integer_,
  SNP_ID = common,
  PIP    = as.numeric(rst$pip),
  stringsAsFactors = FALSE
)

# Overlay signal cluster membership from get_enloc()
enloc <- tryCatch(get_enloc(rst), error = function(e) NULL)
if (!is.null(enloc) && nrow(enloc) > 0) {
  enloc$SNP_ID <- common[as.integer(enloc$SNP_ID)]
  # Update SC_ID for credible set members
  cs_idx <- match(enloc$SNP_ID, all_pips$SNP_ID)
  all_pips$SC_ID[cs_idx[!is.na(cs_idx)]] <- enloc$SC_ID[!is.na(cs_idx)]
}

results[[target_gene]] <- all_pips

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
