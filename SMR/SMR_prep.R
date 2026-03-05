#!/usr/bin/env Rscript
# Select representative SNPs for SMR using mashr results.
# SNP selection priority:
#   1. eGenes: >=1 eQTL with lfsr < 0.05 in either condition
#      -> representative SNP = lowest min(lfsr_trt, lfsr_ctl)
#   2. response eGenes (subset of eGenes): >=1 eQTL with lfsr < 0.05 AND
#      effect differs in direction OR magnitude > fc_threshold
#      -> representative SNP = lowest lfsr among SNPs meeting effect criteria
# Output: SMR input table (with GWAS), SNP classification table

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(annotables)
})

args <- commandArgs(trailingOnly=T)
trait <- ifelse(length(args) > 0, args[1], "CARDIoGRAM_C4D_CAD")

lfsr_threshold <- 0.05
fc_threshold   <- 2
min_effect     <- 0.1  # minimum |posterior mean| to consider directional

mash_base <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/mash_results"
gwas_fn   <- paste0("/nfs/rprdata/julong/sc_multiome/gwas/gwas_prepare/gwas_imputefile/",
                    trait, "_impute_gwas.txt.gz")
outdir    <- paste0("./1_SMR_output/", trait, "/")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)


pairs <- list(
  BPA_T6  = list(trt="BPA_T6",  ctl="H2O_T6",  mash_dir="BPA_H2O_T6"),
  BPA_T24 = list(trt="BPA_T24", ctl="H2O_T24", mash_dir="BPA_H2O_T24"),
  MBP_T6  = list(trt="MBP_T6",  ctl="EtOH_T6",  mash_dir="MBP_EtOH_T6"),
  MBP_T24 = list(trt="MBP_T24", ctl="EtOH_T24", mash_dir="MBP_EtOH_T24")
)

# ── Load GWAS ─────────────────────────────────────────────────────────────────
cat("Loading GWAS...\n")
summ <- fread(gwas_fn, header=TRUE, data.table=FALSE) %>%
  dplyr::select(id_b38, chr, pos, zscore_gwas=zscore, pval_gwas=pval) %>%
  mutate(variant_id = id_b38 %>%
           sub("^chr", "", .) %>%
           sub("_b38$", "", .) %>%
           gsub("_", ":", .))

# Gene symbol lookup via annotables (strips Ensembl version suffix if present)
gene_map <- grch38 %>%
  dplyr::select(ensgene, symbol) %>%
  distinct() %>%
  rename(phenotype_id = ensgene, gene_symbol = symbol)

# ── Helper: load mashr matrix ─────────────────────────────────────────────────
load_mash_mat <- function(path) {
  m <- fread(path, data.table=FALSE)
  rownames(m) <- m[[1]]
  m[, -1, drop=FALSE]
}

# ── Process each pair ─────────────────────────────────────────────────────────
smr_list  <- vector("list", length(pairs))
snp_class_list <- vector("list", length(pairs))

for (pair_name in names(pairs)) {
  trt      <- pairs[[pair_name]]$trt
  ctl      <- pairs[[pair_name]]$ctl
  mash_dir <- file.path(mash_base, pairs[[pair_name]]$mash_dir)
  cat("Processing:", pair_name, "\n")

  # Load mashr outputs
  lfsr <- load_mash_mat(file.path(mash_dir, "mash_lfsr_combined.txt.gz"))
  pm   <- load_mash_mat(file.path(mash_dir, "mash_posterior_mean_combined.txt.gz"))

  # Subset to this pair's conditions
  lfsr_p <- lfsr[, c(trt, ctl), drop=FALSE]
  pm_p   <- pm[,   c(trt, ctl), drop=FALSE]

  # Parse variant_id and phenotype_id from rownames (format: phenotype_id_variant_id)
  # Assumes variant_id contains colons (chr:pos:ref:alt) — split on first underscore
  rn <- rownames(lfsr_p)
  # variant_gene rows: "<gene_id>_<chrom:pos:ref:alt>"
  # Use a regex that splits at the boundary before chr or numeric chrom
  split_idx <- regexpr("_(\\d+:|chr)", rn, perl=TRUE)
  phenotype_id <- substr(rn, 1, split_idx - 1)
  variant_id   <- substr(rn, split_idx + 1, nchar(rn))

  df <- data.frame(
    variant_gene = rn,
    phenotype_id = phenotype_id,
    variant_id   = variant_id,
    lfsr_trt = lfsr_p[[trt]],
    lfsr_ctl = lfsr_p[[ctl]],
    pm_trt   = pm_p[[trt]],
    pm_ctl   = pm_p[[ctl]],
    stringsAsFactors = FALSE
  )

  # ── Classify eQTLs ──────────────────────────────────────────────────────────
  df <- df %>% mutate(
    sig_either = lfsr_trt < lfsr_threshold | lfsr_ctl < lfsr_threshold,
    # direction differs
    opposite_dir = sign(pm_trt) != sign(pm_ctl) &
                   (abs(pm_trt) > min_effect | abs(pm_ctl) > min_effect),
    # magnitude differs by > fc_threshold
    fc = ifelse(abs(pm_ctl) > 0,
                abs(pm_trt / pm_ctl),
                ifelse(abs(pm_trt) > min_effect, Inf, 1)),
    magnitude_differs = fc > fc_threshold | fc < 1/fc_threshold,
    min_effect_either = abs(pm_trt) > min_effect | abs(pm_ctl) > min_effect,
    is_response_eqtl = sig_either & min_effect_either & (opposite_dir | magnitude_differs),
    # smallest lfsr across the pair
    min_lfsr = pmin(lfsr_trt, lfsr_ctl)
  )

  # ── Classify genes ──────────────────────────────────────────────────────────
  # eGenes: any sig_either SNP
  egenes <- df %>%
    filter(sig_either) %>%
    pull(phenotype_id) %>%
    unique()

  # response eGenes: any response eQTL
  regenes <- df %>%
    filter(is_response_eqtl) %>%
    pull(phenotype_id) %>%
    unique()

  cat(sprintf("  eGenes: %d | response eGenes: %d\n",
              length(egenes), length(regenes)))

  # ── Select representative SNP per gene ──────────────────────────────────────
  # For response eGenes: best (lowest lfsr) response eQTL
  rep_regene <- df %>%
    filter(phenotype_id %in% regenes, is_response_eqtl) %>%
    group_by(phenotype_id) %>%
    slice_min(order_by=min_lfsr, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    mutate(gene_class="response_eGene")

  # For eGenes (non-response): best sig eQTL
  rep_egene_only <- df %>%
    filter(phenotype_id %in% setdiff(egenes, regenes), sig_either) %>%
    group_by(phenotype_id) %>%
    slice_min(order_by=min_lfsr, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    mutate(gene_class="eGene")

  rep_snps <- bind_rows(rep_regene, rep_egene_only) %>%
    mutate(
      selected_from = case_when(
        lfsr_trt <= lfsr_ctl ~ "treatment",
        TRUE                 ~ "control"
      ),
      pair = pair_name
    )

  # ── Build SMR output ────────────────────────────────────────────────────────
  smr <- rep_snps %>%
    left_join(summ, by="variant_id") %>%
    filter(!is.na(pval_gwas))

  smr_list[[pair_name]]      <- smr
  snp_class_list[[pair_name]] <- rep_snps
}

smr_all      <- bind_rows(smr_list)
snp_class_all <- bind_rows(snp_class_list)

# ── Attach gene symbols ───────────────────────────────────────────────────────
# Strip Ensembl version suffix (e.g. ENSG00000001234.5 -> ENSG00000001234)
smr_all$phenotype_id       <- sub("\\.\\d+$", "", smr_all$phenotype_id)
snp_class_all$phenotype_id <- sub("\\.\\d+$", "", snp_class_all$phenotype_id)

smr_all       <- left_join(smr_all,       gene_map, by="phenotype_id")
snp_class_all <- left_join(snp_class_all, gene_map, by="phenotype_id")

# ── Write outputs ─────────────────────────────────────────────────────────────
smr_cols <- c("pair", "phenotype_id", "gene_symbol", "gene_class",
              "variant_id", "chr", "pos",
              "lfsr_trt", "lfsr_ctl", "min_lfsr", "selected_from",
              "pm_trt", "pm_ctl",
              "zscore_gwas", "pval_gwas")

class_cols <- c("pair", "phenotype_id", "gene_symbol", "gene_class",
                "variant_id", "selected_from",
                "lfsr_trt", "lfsr_ctl", "min_lfsr",
                "pm_trt", "pm_ctl",
                "is_response_eqtl", "opposite_dir", "magnitude_differs", "fc")

smr_out   <- smr_all[,       intersect(smr_cols,   colnames(smr_all))]
class_out <- snp_class_all[, intersect(class_cols, colnames(snp_class_all))]

smr_fn   <- paste0(outdir, trait, "_mashr_smr_input.txt.gz")
class_fn <- paste0(outdir, trait, "_snp_classification.txt.gz")

write.table(smr_out,   gzfile(smr_fn),   row.names=FALSE, quote=FALSE, sep="\t")
write.table(class_out, gzfile(class_fn), row.names=FALSE, quote=FALSE, sep="\t")

cat("\nOutput written to:\n  ", smr_fn, "\n  ", class_fn, "\n")
cat("SMR input rows:", nrow(smr_out),
    "| SNP classification rows:", nrow(class_out), "\n")
