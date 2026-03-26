#!/usr/bin/env Rscript
# 2_build_eqtl_inputs.R
# Build per-condition TORUS eQTL input files, gene map, and SNP map.
# Uses fwrite() and shell commands throughout for maximum speed.
#
# Outputs (all in OUTDIR):
#   split/{condition}.eQTL.txt.gz  — per-condition z-score file for TORUS
#   zzz_gene.map.gz                — strand-aware TSS per gene
#   zzz_snp.map.gz                 — chromosome and position per SNP

library(data.table)

###############################################################################
# CONFIG
###############################################################################

EQTL_FILE   <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/tensor/combined_nominal/tensorqtl_nominal_SV1-10.txt.gz"
GENCODE_GFF <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz"
OUTDIR      <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_input"

###############################################################################

dir.create(file.path(OUTDIR, "split"), showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 1. Read nominal file once, compute z-score, split by condition
###############################################################################

cat("Reading tensorQTL nominal file...\n")
dat <- fread(EQTL_FILE, header = TRUE)
cat("  Rows:", nrow(dat), "\n")

# Compute z-score and strip gene version in place
dat[, gene   := sub("\\..*", "", phenotype_id)]
dat[, t_stat := slope / slope_se]

cat("  Conditions:", paste(sort(unique(dat$condition)), collapse = ", "), "\n\n")

for (cond in sort(unique(dat$condition))) {
  cat("  Writing:", cond)
  out <- dat[condition == cond, .(SNP = variant_id, gene, beta = slope,
                                   `t-stat` = t_stat, `p-value` = pval_nominal)]
  fwrite(out, file.path(OUTDIR, "split", paste0(cond, ".eQTL.txt.gz")),
         sep = " ", compress = "gzip")
  cat(" ->", nrow(out), "rows\n")
}

###############################################################################
# 2. SNP map: SNP chromosome position
###############################################################################

cat("\nBuilding SNP map...\n")

snp_map <- unique(dat[, .(SNP = variant_id)])
snp_map[, c("chr", "pos", "ref", "alt") := tstrsplit(SNP, ":", fixed = TRUE)]
snp_map <- snp_map[, .(SNP, chr, pos)]

cat("  Unique SNPs:", nrow(snp_map), "\n")
fwrite(snp_map, file.path(OUTDIR, "zzz_snp.map.gz"),
       sep = " ", col.names = FALSE, compress = "gzip")
cat("  Written: zzz_snp.map.gz\n")

###############################################################################
# 3. Gene map: strand-aware TSS from Gencode v31
###############################################################################

cat("\nBuilding gene map from Gencode v31...\n")

gene_list <- unique(dat[, .(phenotype_id, gene_base = gene)])
cat("  Genes in eQTL data:", nrow(gene_list), "\n")

cat("  Reading Gencode GFF3...\n")
anno <- fread(
  cmd = paste("zcat", GENCODE_GFF, "| grep -v '^#' | awk '$3==\"gene\"'"),
  header = FALSE, sep = "\t",
  col.names = c("chr", "source", "feature", "start", "end",
                "score", "strand", "frame", "attributes")
)

anno[, gene_id_full := sub("ID=([^;]+).*", "\\1", attributes)]
anno[, gene_base    := sub("\\..*", "", gene_id_full)]
anno[, chr_clean    := sub("chr", "", chr)]
anno[, tss          := ifelse(strand == "+", start, end)]
anno <- anno[!grepl("PAR_Y", gene_id_full)]

gene_map <- merge(gene_list, anno[, .(gene_base, chr_clean, tss)],
                  by = "gene_base", all.x = FALSE)
gene_map <- gene_map[!is.na(tss), .(gene = phenotype_id, chr = chr_clean,
                                     start = tss, start2 = tss)]

cat("  Genes mapped:", nrow(gene_map), "\n")
cat("  Unmapped:", nrow(gene_list) - nrow(gene_map), "\n")

fwrite(gene_map, file.path(OUTDIR, "zzz_gene.map.gz"),
       sep = "\t", col.names = FALSE, compress = "gzip")
cat("  Written: zzz_gene.map.gz\n")

cat("\nDone.\n")
