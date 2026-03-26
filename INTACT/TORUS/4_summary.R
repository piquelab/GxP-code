#!/usr/bin/env Rscript
# 4_summary.R
# Parse TORUS enrichment results and produce forest plots.
#
# Outputs:
#   enrichment_summary.txt           — full table
#   enrichment_by_condition.png      — conditions on Y, annotations as colors
#   enrichment_by_annotation.png     — annotations on Y, conditions as colors
#   enrichment_combined.png          — both panels stacked

library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)

###############################################################################
# CONFIG
###############################################################################

OUTDIR     <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_output"
CONDITIONS <- c("BPA_T6", "BPA_T24", "EtOH_T6", "EtOH_T24",
                "H2O_T6", "H2O_T24", "MBP_T6",  "MBP_T24")
COND_ORDER <- c("H2O_T6", "H2O_T24", "EtOH_T6", "EtOH_T24",
                "BPA_T6", "BPA_T24", "MBP_T6",  "MBP_T24")

###############################################################################
# Parse .est files
###############################################################################

cat("Parsing .est files...\n")
results <- lapply(CONDITIONS, function(cond) {
  fn <- file.path(OUTDIR, paste0(cond, ".est"))
  if (!file.exists(fn)) { cat("  Missing:", fn, "\n"); return(NULL) }
  x <- read.table(fn, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("annotation", "estimate", "CI_lower", "CI_upper"))
  x$condition <- cond
  x
})
est <- rbindlist(Filter(Negate(is.null), results))

write.table(est, file.path(OUTDIR, "enrichment_summary.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
cat("Written: enrichment_summary.txt\n\n")

# Extract bin number safely
est$bin_num <- suppressWarnings(
  as.integer(sub(".*dtss\\.([0-9]+).*", "\\1", est$annotation))
)

# Filter: open chromatin + dtss bins 4-17 only
est <- est[grepl("open_chromatin", est$annotation) |
           (!is.na(est$bin_num) & est$bin_num >= 4 & est$bin_num <= 17), ]

# Clean labels
est$label <- ifelse(
  grepl("open_chromatin", est$annotation),
  "Open chromatin",
  gsub("dtss\\.[0-9]+::", "", est$annotation)
)

# Factor levels: open chromatin first, then TSS bins proximal to distal
tss_labels <- unique(est$label[est$label != "Open chromatin"])
tss_order  <- tss_labels[order(est$bin_num[match(tss_labels, est$label)])]
annot_order <- c("Open chromatin", tss_order)
est$label     <- factor(est$label, levels = annot_order)
est$condition <- factor(est$condition, levels = COND_ORDER)

# Numeric y for manual dodging
est$y_cond  <- as.integer(est$condition)
est$y_annot <- as.integer(est$label)

n_annot <- nlevels(est$label)
n_cond  <- nlevels(est$condition)

# Dodge offset per group
dodge_range <- 0.7
make_offsets <- function(n) seq(-dodge_range/2, dodge_range/2, length.out = n)

###############################################################################
# Color palettes
###############################################################################

n_tss   <- length(tss_order)
tss_pal <- colorRampPalette(c("#B5D4F4", "#042C53"))(n_tss)
names(tss_pal) <- tss_order
annot_pal <- c("Open chromatin" = "black", tss_pal)

cond_pal <- c(
  "H2O_T6"  = "#B4B2A9", "H2O_T24"  = "#5F5E5A",
  "EtOH_T6" = "#888780", "EtOH_T24" = "#2C2C2A",
  "BPA_T6"  = "#85B7EB", "BPA_T24"  = "#185FA5",
  "MBP_T6"  = "#F0997B", "MBP_T24"  = "#993C1D"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.y       = element_blank(),
    plot.title         = element_text(size = 11, face = "bold")
  )

###############################################################################
# Plot A: conditions on Y, annotations as colors
###############################################################################

# Add dodge offset per annotation within each condition row
annot_offsets <- make_offsets(n_annot)
names(annot_offsets) <- levels(est$label)
est_a <- est %>%
  mutate(y_dodge = y_cond + annot_offsets[as.character(label)])

p_a <- ggplot(est_a, aes(color = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, y = y_dodge),
                 height = 0, linewidth = 0.45, alpha = 0.7) +
  geom_point(aes(x = estimate, y = y_dodge), size = 1.8) +
  scale_color_manual("Annotation", values = annot_pal,
                     guide = guide_legend(ncol = 2,
                                          override.aes = list(size = 2.5))) +
  scale_x_continuous("Log odds ratio (95% CI)") +
  scale_y_continuous(breaks = seq_len(n_cond),
                     labels = levels(est$condition)) +
  ggtitle("TORUS enrichment by condition") +
  base_theme +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 7.5),
        legend.title    = element_text(size = 9),
        legend.key.size = unit(0.75, "lines"))

ggsave(file.path(OUTDIR, "enrichment_by_condition.png"), p_a,
       width = 7.0, height = 5.5, dpi = 150)
cat("Written: enrichment_by_condition.png\n")

###############################################################################
# Plot B: annotations on Y, conditions as colors
###############################################################################

cond_offsets <- make_offsets(n_cond)
names(cond_offsets) <- levels(est$condition)
est_b <- est %>%
  mutate(y_dodge = y_annot + cond_offsets[as.character(condition)])

p_b <- ggplot(est_b, aes(color = condition)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, y = y_dodge),
                 height = 0, linewidth = 0.45, alpha = 0.7) +
  geom_point(aes(x = estimate, y = y_dodge), size = 1.8) +
  scale_color_manual("Condition", values = cond_pal,
                     guide = guide_legend(ncol = 2,
                                          override.aes = list(size = 2.5))) +
  scale_x_continuous("Log odds ratio (95% CI)") +
  scale_y_continuous(breaks = seq_len(n_annot),
                     labels = levels(est$label)) +
  ggtitle("TORUS enrichment by annotation") +
  base_theme +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9),
        legend.title    = element_text(size = 9),
        legend.key.size = unit(0.75, "lines"))

ggsave(file.path(OUTDIR, "enrichment_by_annotation.png"), p_b,
       width = 7.0, height = 6.5, dpi = 150)
cat("Written: enrichment_by_annotation.png\n")

###############################################################################
# Combined: stacked
###############################################################################

p_combined <- plot_grid(p_a, p_b, ncol = 1, align = "v",
                        labels = c("A", "B"), label_size = 12)
ggsave(file.path(OUTDIR, "enrichment_combined.png"), p_combined,
       width = 7.5, height = 12.0, dpi = 150)
cat("Written: enrichment_combined.png\n")

cat("\nDone.\n")
