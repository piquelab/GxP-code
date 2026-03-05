#!/usr/bin/env Rscript
# Summarize SMR results
# Ouputs:
#   1. Summary table with significant genes per pair at FDR < 0.10
#   2. P-value distributions 
#   3. Manhattan plots, annotated with FDR < 0.01

library(tidyverse)
library(data.table)
library(qvalue)
library(annotables)
library(openxlsx)
library(ggrepel)
rm(list=ls())

trait <- "CARDIoGRAM_C4D_CAD"

pairs <- list(
  BPA_T6  = c(treatment="BPA_T6",  control="H2O_T6"),
  BPA_T24 = c(treatment="BPA_T24", control="H2O_T24"),
  MBP_T6  = c(treatment="MBP_T6",  control="EtOH_T6"),
  MBP_T24 = c(treatment="MBP_T24", control="EtOH_T24")
)

# Significance thresholds
FDR_SIG   <- 0.10  # used for calling significant genes throughout
FDR_LABEL <- 0.01  # used only for Manhattan labels

pair_cols <- c(BPA_T6="#4dac26", BPA_T24="#b8e186", MBP_T6="#d01c8b", MBP_T24="#f1b6da")

grch38_unq <- grch38 %>%
  filter(chr %in% as.character(1:22), grepl("protein", biotype)) %>%
  distinct(ensgene, .keep_all=TRUE) %>%
  dplyr::select(gene=ensgene, symbol)

outdir <- "./2_summary_twas/"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

###############################################################################
### Part 1: FDR calculation and per-pair output
###############################################################################

smr <- fread(paste0("./1_SMR_output/", trait, "/", trait, "_mashr_smr_input.txt.gz"),
             header=TRUE, data.table=FALSE)

for (pair_name in names(pairs)) {
  outdir2 <- paste0(outdir, pair_name, "/")
  dir.create(outdir2, showWarnings=FALSE, recursive=TRUE)

  res <- smr %>%
    filter(pair == pair_name) %>%
    mutate(FDR = qvalue(pval_gwas)$qvalues) %>%
    left_join(grch38_unq, by=c("phenotype_id"="gene")) %>%
    mutate(symbol = ifelse(is.na(symbol), "", symbol)) %>%
    arrange(pval_gwas)

  write.table(res, paste0(outdir2, trait, "_", pair_name, "_twas.txt"),
              row.names=FALSE, col.names=TRUE, quote=FALSE, na="")
  cat(pair_name, "—", sum(res$FDR < FDR_SIG), "genes at FDR <", FDR_SIG, "\n")
}

###############################################################################
### Part 2: Summary table — significant genes per pair at FDR < 0.10
###############################################################################

summ <- map_dfr(names(pairs), function(pair_name) {
  x <- fread(paste0(outdir, pair_name, "/", trait, "_", pair_name, "_twas.txt"),
             header=TRUE, data.table=FALSE, fill=TRUE)
  data.frame(
    pair         = pair_name,
    n_sig_FDR10  = x %>% filter(FDR < FDR_SIG) %>% pull(phenotype_id) %>% n_distinct()
  )
})
write.xlsx(summ, paste0(outdir, "1_summary_ngene_FDR0.10.xlsx"), overwrite=TRUE)

###############################################################################
### Part 3: Proportion of genes with p < 0.01 per pair
###############################################################################

summ_prop <- map_dfr(names(pairs), function(pair_name) {
  x <- smr %>% filter(pair == pair_name)
  data.frame(pair=pair_name, prop_p01=round(mean(x$pval_gwas < 0.01), digits=3))
})
write.xlsx(summ_prop, paste0(outdir, "1.2_summary_prop_sig0.01.xlsx"), overwrite=TRUE)

###############################################################################
### Part 4: P-value histograms with pi0 annotation
###############################################################################

outdir_hist <- paste0(outdir, "hist/")
dir.create(outdir_hist, showWarnings=FALSE, recursive=TRUE)

for (pair_name in names(pairs)) {
  col <- pair_cols[pair_name]

  res <- fread(paste0(outdir, pair_name, "/", trait, "_", pair_name, "_twas.txt"),
               header=TRUE, data.table=FALSE, fill=TRUE) %>%
    filter(!is.na(pval_gwas))

  fdr_res <- qvalue(res$pval_gwas)
  pi0     <- round(fdr_res$pi0, digits=3)
  n_bins  <- 40
  null_y  <- (nrow(res) * pi0) / n_bins

  p <- ggplot(res, aes(x=pval_gwas)) +
    geom_histogram(bins=n_bins, fill=col, color=NA, alpha=0.75) +
    geom_hline(yintercept=null_y, color=col, linetype="dashed", linewidth=0.7) +
    annotate("text", x=0.75, y=null_y * 1.15,
             label=bquote(pi[0]==.(pi0)), size=4.5, color="grey20") +
    scale_x_continuous("TWAS p-value", expand=c(0, 0)) +
    scale_y_continuous("Number of genes", expand=expansion(mult=c(0, 0.2))) +
    ggtitle(pair_name) +
    theme_bw(base_size=12) +
    theme(
      plot.title     = element_text(hjust=0.5, face="bold"),
      panel.grid     = element_blank(),
      axis.line      = element_line(color="grey30"),
      plot.margin    = margin(10, 15, 10, 10)
    )

  ggsave(paste0(outdir_hist, "Figure1_", trait, "_", pair_name, "_pval_hist.png"),
         p, width=5, height=3.5, dpi=150)
  cat("Histogram saved:", pair_name, "\n")
}

###############################################################################
### Part 5: Manhattan plots
# Significant threshold (dashed line): FDR < FDR_SIG (0.10)
# Labels: FDR < FDR_LABEL (0.01) only
###############################################################################

for (pair_name in names(pairs)) {
  col <- pair_cols[pair_name]

  res <- fread(paste0(outdir, pair_name, "/", trait, "_", pair_name, "_twas.txt"),
               header=TRUE, data.table=FALSE, fill=TRUE) %>%
    filter(!is.na(pval_gwas), !is.na(chr), !is.na(pos)) %>%
    mutate(
      chr    = as.integer(chr),
      log10p = -log10(pval_gwas),
      gr2    = factor(chr %% 2),
      is_sig = FDR < FDR_SIG
    ) %>%
    arrange(chr, pos)

  # Cumulative genomic position
  chrom_offsets <- res %>%
    group_by(chr) %>%
    summarise(max_pos=max(pos), .groups="drop") %>%
    arrange(chr) %>%
    mutate(offset = lag(cumsum(max_pos), default=0))

  res <- res %>%
    left_join(chrom_offsets[, c("chr", "offset")], by="chr") %>%
    mutate(BPcum = pos + offset)

  axis_set <- res %>%
    group_by(chr) %>%
    summarise(midpoint=(max(BPcum) + min(BPcum)) / 2, .groups="drop")

  # Threshold line at FDR_SIG
  sig_line <- res %>% filter(is_sig) %>% pull(log10p) %>% min(na.rm=TRUE)
  sig_line <- if (is.finite(sig_line)) sig_line else NA

  # Labels at FDR_LABEL
  label_df <- res %>% filter(FDR < FDR_LABEL)

  p <- ggplot(res, aes(x=BPcum, y=log10p)) +
    geom_point(aes(color=gr2, size=log10p), alpha=0.8, shape=16) +
    { if (!is.na(sig_line))
        geom_hline(yintercept=sig_line, color=col, linetype="dashed", linewidth=0.6) } +
    { if (nrow(label_df) > 0)
        geom_text_repel(data=label_df, aes(label=symbol),
                        fontface="italic", size=2.8,
                        box.padding=0.3, point.padding=0.2,
                        max.overlaps=30, segment.color="grey50",
                        segment.size=0.3, min.segment.length=0.2) } +
    scale_x_continuous("Chromosome", labels=axis_set$chr, breaks=axis_set$midpoint,
                       expand=expansion(mult=0.01)) +
    scale_y_continuous(bquote(-log[10](italic(p))),
                       limits=c(0, max(res$log10p, na.rm=TRUE) + 3),
                       expand=expansion(mult=c(0, 0.05))) +
    scale_color_manual(values=c("0"="#2171b5", "1"="#6baed6")) +
    scale_size_continuous(range=c(0.3, 2.0)) +
    ggtitle(pair_name) +
    theme_bw(base_size=12) +
    theme(
      legend.position   = "none",
      panel.grid        = element_blank(),
      axis.text.x       = element_text(size=9, angle=0, vjust=0.5),
      axis.text.y       = element_text(size=10),
      axis.title        = element_text(size=11),
      plot.title        = element_text(hjust=0.5, face="bold"),
      plot.margin       = margin(10, 15, 10, 10)
    )

  ggsave(paste0(outdir, "Figure2_", trait, "_", pair_name, "_manhattan.png"),
         p, width=9, height=4.5, dpi=150)
  cat("Manhattan saved:", pair_name, "\n")
}
