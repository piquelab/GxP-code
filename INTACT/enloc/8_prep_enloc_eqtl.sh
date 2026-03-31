#!/usr/bin/env Rscript
# 8_prep_enloc_eqtl.R
# Convert DAP-S PIP output to fastENLOC eQTL annotation format.
# One file per condition.
#
# fastENLOC v3 eQTL format (no header, 6 tab-separated columns):
#   CHROM  POS  SNP_ID  REF  ALT  ANNOTATION
#
# ANNOTATION format per variant:
#   gene:signal_cluster@=variant_PIP[cluster_PIP:n_variants_in_cluster]
#   Multiple annotations separated by ;
#
# Example:
#   chr1 54421 chr1_54421_A_G_b38 A G ENSG00000227232:3@=1.05e-04[2.66e-01:2]

library(data.table)

###############################################################################
# CONFIG
###############################################################################

DAPSDIR <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/daps_output"
OUTDIR  <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/enloc_input"
CONDITIONS <- c("BPA_T6", "BPA_T24", "MBP_T6", "MBP_T24")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

###############################################################################

for (cond in CONDITIONS) {
  cat("Processing:", cond, "\n")

  pip <- fread(file.path(DAPSDIR, paste0(cond, "_daps_pips.txt.gz")),
               header = TRUE)
  setnames(pip, c("gene", "SC_ID", "SNP_ID", "PIP"))

  # Drop variants not in any signal cluster
  pip <- pip[!is.na(SC_ID)]
  cat("  Variants in signal clusters:", nrow(pip), "\n")

  # Parse SNP coordinates
  pip[, c("chr_num", "pos", "ref", "alt") := tstrsplit(SNP_ID, ":", fixed = TRUE)]
  pip[, pos := as.integer(pos)]

  # Convert SNP ID to chrN_pos_ref_alt_b38 format
  pip[, snp_b38 := paste0("chr", chr_num, "_", pos, "_", ref, "_", alt, "_b38")]
  pip[, chrom   := paste0("chr", chr_num)]

  # Compute cluster-level summary stats
  pip[, n_var      := .N,        by = .(gene, SC_ID)]
  pip[, cluster_pip := max(PIP), by = .(gene, SC_ID)]

  # Build annotation string per variant:
  # gene:SC_ID@=variant_PIP[cluster_PIP:n_variants]
  pip[, annot := paste0(
    gene, ":", SC_ID,
    "@=", formatC(PIP,         format = "e", digits = 5),
    "[",  formatC(cluster_pip, format = "e", digits = 3),
    ":",  n_var, "]"
  )]

  # A SNP can appear in multiple gene:cluster annotations — collapse with ;
  snp_annot <- pip[, .(
    chrom = chrom[1],
    pos   = pos[1],
    ref   = ref[1],
    alt   = alt[1],
    annot = paste(annot, collapse = ";")
  ), by = snp_b38]

  setorder(snp_annot, chrom, pos)

  # Write — no header, tab-separated
  out <- snp_annot[, .(chrom, pos, snp_b38, ref, alt, annot)]
  outfile <- file.path(OUTDIR, paste0(cond, "_enloc.eqtl.vcf.gz"))
  fwrite(out, outfile, sep = "\t", col.names = FALSE, compress = "gzip",
         quote = FALSE)

  cat("  Written:", outfile, "\n")
  cat("  Unique SNPs:", nrow(out),
      "| Genes:", uniqueN(pip$gene),
      "| Clusters:", uniqueN(paste(pip$gene, pip$SC_ID)), "\n\n")
}

cat("Done.\n")
