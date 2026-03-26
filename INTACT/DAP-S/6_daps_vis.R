#!/usr/bin/env Rscript
# 6_daps_vis.R
# Visualize DAP-S fine-mapping results across 8 conditions.
#
# Plots produced:
#   1. pip_distribution.png       — PIP histogram across all conditions
#   2. credible_set_size.png      — signal cluster size distribution
#   3. high_pip_counts.png        — variants with PIP > 0.5 per condition
#   4. example_locus.png          — per-SNP PIPs for highest-confidence gene

library(data.table)
library(ggplot2)
library(cowplot)

###############################################################################
# CONFIG
###############################################################################

DAPSDIR <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/daps_output"
OUTDIR  <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/daps_output/figures"
CONDITIONS <- c("BPA_T6", "BPA_T24", "EtOH_T6", "EtOH_T24",
                "H2O_T6", "H2O_T24", "MBP_T6",  "MBP_T24")
COND_ORDER <- c("H2O_T6", "H2O_T24", "EtOH_T6", "EtOH_T24",
                "BPA_T6", "BPA_T24", "MBP_T6",  "MBP_T24")
PIP_THRESH <- 0.5

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Load all PIP data
###############################################################################

cat("Loading PIP data...\n")
dat <- rbindlist(lapply(CONDITIONS, function(cond) {
  f <- file.path(DAPSDIR, paste0(cond, "_daps_pips.txt.gz"))
  x <- fread(f, header = TRUE)
  x$condition <- cond
  x
}))
setnames(dat, c("gene", "SC_ID", "SNP_ID", "PIP", "condition"))
dat$condition <- factor(dat$condition, levels = COND_ORDER)
dat$group <- ifelse(dat$condition %in% c("H2O_T6","H2O_T24","EtOH_T6","EtOH_T24"),
                    "Control", "Treatment")

cat("  Total PIP rows:", nrow(dat), "\n")
cat("  Unique genes:", uniqueN(dat$gene), "\n")
cat("  Unique variants:", uniqueN(dat$SNP_ID), "\n\n")

base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        plot.title        = element_text(size = 11, face = "bold"),
        plot.subtitle     = element_text(size = 9, color = "grey40"),
        strip.background  = element_rect(fill = "grey95"),
        strip.text        = element_text(size = 9))

###############################################################################
# Plot 1: PIP distribution — faceted by condition
###############################################################################

cat("Plot 1: PIP distribution...\n")

p1 <- ggplot(dat, aes(x = PIP)) +
  geom_histogram(bins = 50, fill = "grey40", color = NA) +
  geom_vline(xintercept = PIP_THRESH, linetype = "dashed",
             color = "#E24B4A", linewidth = 0.5) +
  facet_wrap(~condition, nrow = 2) +
  scale_x_continuous("Posterior inclusion probability (PIP)",
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous("Number of variants") +
  annotate("text", x = 0.52, y = Inf, label = paste("PIP >", PIP_THRESH),
           hjust = 0, vjust = 1.5, size = 2.8, color = "#E24B4A") +
  ggtitle("DAP-S fine-mapping: PIP distribution",
          subtitle = "Most variants have low PIP; high-confidence causal variants cluster near 1") +
  base_theme

ggsave(file.path(OUTDIR, "pip_distribution.png"), p1,
       width = 8, height = 5, dpi = 150)
cat("  Written: pip_distribution.png\n")

###############################################################################
# Plot 2: Credible set (signal cluster) size distribution
###############################################################################

cat("Plot 2: Credible set sizes...\n")

# Count variants per signal cluster per gene per condition
cs_size <- dat[!is.na(SC_ID),
               .(size = .N),
               by = .(condition, gene, SC_ID)]

# Summary stats for subtitle
med_size <- median(cs_size$size)
pct_small <- round(100 * mean(cs_size$size <= 10))

p2 <- ggplot(cs_size, aes(x = pmin(size, 50))) +
  geom_histogram(bins = 40, fill = "grey40", color = NA) +
  facet_wrap(~condition, nrow = 2) +
  scale_x_continuous("Signal cluster size (variants; capped at 50)",
                     breaks = c(1, 10, 20, 30, 40, 50)) +
  scale_y_continuous("Number of signal clusters") +
  ggtitle("DAP-S fine-mapping: signal cluster sizes",
          subtitle = sprintf("Median size = %d variants; %d%% of clusters have ≤10 variants",
                             med_size, pct_small)) +
  base_theme

ggsave(file.path(OUTDIR, "credible_set_size.png"), p2,
       width = 8, height = 5, dpi = 150)
cat("  Written: credible_set_size.png\n")

###############################################################################
# Plot 3: High-confidence variants (PIP > threshold) per condition
###############################################################################

cat("Plot 3: High-PIP variant counts...\n")

high_pip <- dat[PIP > PIP_THRESH,
                .(n_variants = uniqueN(SNP_ID),
                  n_genes    = uniqueN(gene)),
                by = condition]
high_pip$condition <- factor(high_pip$condition, levels = COND_ORDER)

# Two-panel: variants and genes
p3a <- ggplot(high_pip, aes(x = condition, y = n_variants, fill = condition)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_variants), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(
    "H2O_T6" = "#B4B2A9", "H2O_T24" = "#5F5E5A",
    "EtOH_T6" = "#888780", "EtOH_T24" = "#2C2C2A",
    "BPA_T6" = "#85B7EB", "BPA_T24" = "#185FA5",
    "MBP_T6" = "#F0997B", "MBP_T24" = "#993C1D"
  )) +
  scale_y_continuous("Unique variants") +
  ggtitle(paste0("Variants with PIP > ", PIP_THRESH)) +
  base_theme +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

p3b <- ggplot(high_pip, aes(x = condition, y = n_genes, fill = condition)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_genes), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(
    "H2O_T6" = "#B4B2A9", "H2O_T24" = "#5F5E5A",
    "EtOH_T6" = "#888780", "EtOH_T24" = "#2C2C2A",
    "BPA_T6" = "#85B7EB", "BPA_T24" = "#185FA5",
    "MBP_T6" = "#F0997B", "MBP_T24" = "#993C1D"
  )) +
  scale_y_continuous("Unique genes (eGenes)") +
  ggtitle(paste0("Genes with ≥1 variant PIP > ", PIP_THRESH)) +
  base_theme +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- plot_grid(p3a, p3b, nrow = 1, labels = c("A", "B"), label_size = 11)
ggsave(file.path(OUTDIR, "high_pip_counts.png"), p3,
       width = 8, height = 4.5, dpi = 150)
cat("  Written: high_pip_counts.png\n")

###############################################################################
# Plot 4: Example locus — gene with highest max PIP, BPA_T6
###############################################################################

cat("Plot 4: Example locus...\n")

# Pick the gene with the highest single-variant PIP in BPA_T6
ex_cond <- "BPA_T6"
ex_gene <- dat[condition == ex_cond][which.max(PIP), gene]
ex_dat  <- dat[condition == ex_cond & gene == ex_gene]

# Parse chromosome and position from SNP_ID (chr:pos:ref:alt)
ex_dat[, c("chr", "pos", "ref", "alt") := tstrsplit(SNP_ID, ":", fixed = TRUE)]
ex_dat[, pos := as.integer(pos)]
ex_dat[, in_cs := !is.na(SC_ID)]

cat("  Example gene:", ex_gene, "| max PIP:", round(max(ex_dat$PIP), 3), "\n")

p4 <- ggplot(ex_dat, aes(x = pos / 1e6, y = PIP, color = in_cs)) +
  geom_segment(aes(xend = pos / 1e6, yend = 0),
               linewidth = 0.4, alpha = 0.7) +
  geom_point(size = 2) +
  scale_color_manual("In signal cluster",
                     values = c("TRUE" = "#185FA5", "FALSE" = "grey60")) +
  scale_x_continuous("Genomic position (Mb)") +
  scale_y_continuous("Posterior inclusion probability (PIP)",
                     limits = c(0, 1.05)) +
  ggtitle(paste0("Example locus: ", ex_gene, " (", ex_cond, ")"),
          subtitle = paste0("Variants colored by signal cluster membership; ",
                            "max PIP = ", round(max(ex_dat$PIP), 3))) +
  base_theme +
  theme(legend.position = "bottom")

ggsave(file.path(OUTDIR, "example_locus.png"), p4,
       width = 7, height = 4, dpi = 150)
cat("  Written: example_locus.png\n")

###############################################################################
# Summary stats to console
###############################################################################

cat("\n--- Summary ---\n")
cat("Total variants tested:       ", nrow(dat), "\n")
cat("Variants with PIP > 0.5:     ", nrow(dat[PIP > 0.5]), "\n")
cat("Variants with PIP > 0.9:     ", nrow(dat[PIP > 0.9]), "\n")
cat("Genes with any PIP > 0.5:    ", uniqueN(dat[PIP > 0.5]$gene), "\n")
cat("Median signal cluster size:  ", median(cs_size$size), "\n")
cat("Done.\n")
