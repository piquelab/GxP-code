#!/usr/bin/env Rscript
# 1_build_annotation.R
# Build zzz_torus.annot.gz for TORUS from either CENTIPEDE footprint BED files
# or a single open chromatin BED file.
#
# CENTIPEDE mode: reads SNP IDs directly from column 9 of each motif BED file.
# Open chromatin mode: uses bedtools intersect to find SNPs overlapping regions.
#
# Output: INTACT/torus_input/zzz_torus.annot.gz
#   Columns: SNP  peaking_d  {motif1}_d  {motif2}_d  ...
#   All annotation columns are binary (0/1), discrete (_d) format for TORUS.

library(tidyverse)
library(data.table)

###############################################################################
# CONFIG — edit this block only
###############################################################################

MODE <- "open_chromatin"   # "centipede" or "open_chromatin"

# CENTIPEDE: directory of per-motif *.bed.gz files
CENTIPEDE_DIR <- "/rs/rs_grp_gxp/hp8265_ATAC_output/Footprinting/CENTIPEDE/Merged_AllBatch_All_Strict_Motifs_SNP_filtered"

# Open chromatin: single BED file (chr start end) — hg38, col 1-3 used
OPEN_CHROMATIN_BED <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/TORUS/open_chromatin_hg38_clean.bed"

# eQTL nominal file — used to get the universe of SNPs to annotate
EQTL_FILE <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/tensor/combined_nominal/tensorqtl_nominal_SV1-10.txt.gz"

# Minimum number of SNPs a motif must annotate to be kept (filters rare annotations)
# 67 motifs total; MA0605.3 (264) and MA0656.2 (266) will be dropped at this threshold
MIN_SNP <- 500

# Output directory
OUTDIR <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_input"

###############################################################################

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("Reading SNP universe from tensorQTL nominal file...\n")
# Read variant_id column directly — fread handles gz natively and reads the full file
snp_universe <- fread(
  EQTL_FILE,
  select = "variant_id",
  header = TRUE,
  data.table = FALSE
)$variant_id |> unique()
cat("  Unique SNPs in eQTL data:", length(snp_universe), "\n")
 
###############################################################################
# Mode: CENTIPEDE
###############################################################################
 
if (MODE == "centipede") {
 
  cat("\nMode: CENTIPEDE\n")
  bed_files <- list.files(CENTIPEDE_DIR, pattern = "\\.bed\\.gz$", full.names = TRUE)
  cat("  Motif files found:", length(bed_files), "\n\n")
 
  # Initialize annotation matrix — one row per SNP in universe
  annot <- data.frame(SNP = snp_universe, stringsAsFactors = FALSE)
 
  for (f in bed_files) {
    motif <- sub("\\.bed\\.gz$", "", basename(f))
    cat("  Processing:", motif, "")
 
    # Column 9 is the SNP ID (chr:pos:ref:alt) — no header in these files
    snps_in_motif <- tryCatch({
      fread(
        cmd = paste("zcat", f, "| awk '{print $9}'"),
        header = FALSE, data.table = FALSE
      )$V1 |> unique()
    }, error = function(e) {
      cat("  [ERROR reading file, skipping]\n")
      return(character(0))
    })
 
    # Intersect with SNP universe
    snps_annotated <- intersect(snps_in_motif, snp_universe)
    n <- length(snps_annotated)
    cat("->", n, "SNPs\n")
 
    if (n < MIN_SNP) {
      cat("    Skipping (below MIN_SNP =", MIN_SNP, ")\n")
      next
    }
 
    col_name <- paste0(motif, "_d")
    annot[[col_name]] <- as.integer(annot$SNP %in% snps_annotated)
  }
 
###############################################################################
# Mode: Open chromatin
###############################################################################
 
} else if (MODE == "open_chromatin") {
 
  cat("\nMode: Open chromatin\n")
  cat("  BED file:", OPEN_CHROMATIN_BED, "\n")
 
  # Convert SNP universe to point BED (0-based start, 1-based end)
  # variant_id format: chr:pos:ref:alt (pos is 1-based)
  snp_bed_file <- file.path(OUTDIR, "tmp_snps.bed")
 
  cat("  Writing SNP BED file...\n")
  snp_df <- data.frame(SNP = snp_universe) |>
    separate(SNP, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) |>
    mutate(
      chr = paste0("chr", chr),   # add chr prefix to match chromatin BED
      start = as.integer(pos) - 1L,
      end = as.integer(pos)
    ) |>
    select(chr, start, end, SNP)
 
 options(scipen = 999)  # prevent scientific notation in BED coordinates
  write.table(snp_df, file = snp_bed_file, quote = FALSE,
              row.names = FALSE, col.names = FALSE, sep = "\t")
 
  # Run bedtools intersect
  intersect_file <- file.path(OUTDIR, "tmp_chromatin_intersect.bed")
  cmd <- paste(
    "bedtools intersect -wa -a", snp_bed_file,
    "-b", OPEN_CHROMATIN_BED, ">", intersect_file
  )
  cat("  Running bedtools intersect...\n")
  system(cmd)
 
  # Read intersected SNPs
  hits <- fread(intersect_file, header = FALSE, data.table = FALSE)
  snps_annotated <- unique(hits$V4)
  cat("  SNPs in open chromatin:", length(snps_annotated), "\n")
 
  # Build annotation matrix
  annot <- data.frame(SNP = snp_universe, stringsAsFactors = FALSE)
  annot[["open_chromatin_d"]] <- as.integer(annot$SNP %in% snps_annotated)
 
  # Clean up temp files
  file.remove(snp_bed_file, intersect_file)
 
} else {
  stop("MODE must be 'centipede' or 'open_chromatin'")
}
 
###############################################################################
# Add peaking_d (union of all annotations) and write output
###############################################################################
 
annot_cols <- setdiff(names(annot), "SNP")
 
if (length(annot_cols) == 0) {
  stop("No annotations passed the MIN_SNP threshold. Lower MIN_SNP or check input files.")
}
 
cat("\nBuilding peaking_d (union of all annotations)...\n")
annot[["peaking_d"]] <- as.integer(rowSums(annot[, annot_cols, drop = FALSE]) > 0)
 
# Reorder: SNP, peaking_d, then individual annotations
annot <- annot[, c("SNP", "peaking_d", annot_cols)]
 
cat("Annotation summary:\n")
cat("  SNPs total:", nrow(annot), "\n")
cat("  SNPs with any annotation (peaking_d=1):", sum(annot$peaking_d), "\n")
cat("  Annotation columns:", length(annot_cols), "\n")
cat("  Column names (first 10):", paste(head(annot_cols, 10), collapse = ", "), "\n")
 
outfile <- file.path(OUTDIR, "zzz_torus_oc.annot.gz")
cat("\nWriting:", outfile, "\n")
fwrite(annot, outfile, quote = FALSE, sep = " ")
cat("Done.\n")
