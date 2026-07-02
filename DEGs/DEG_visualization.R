#!/usr/bin/env Rscript
################################################################################
# Output: figures/ — barplot PDF, volcano+MA PDF
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

################################################################################
# CONFIG
################################################################################

BASE_DIR  <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates"
DEG_FILE  <- file.path(BASE_DIR, "DEG/04202026_AllDEGs_FDR10.txt")
GENE_INFO <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/starter_files/geneInfo.txt"
OUT_DIR   <- file.path(BASE_DIR, "DEG", "figures")
dir.create(OUT_DIR, showWarnings = FALSE)

# Full per-condition results (all tested genes)
ALLGENE_FILES <- list(
  "MBP_T6"  = file.path(BASE_DIR, "DEG/04202026_MBP_T6.txt"),
  "MBP_T24" = file.path(BASE_DIR, "DEG/04202026_MBP_T24.txt")
)

FDR_THRESH <- 0.10
FC_THRESH  <- 0       # for bar plot DEG counting
FC_LINE    <- 0.1     # vertical lines on volcano
N_LABEL    <- 15
TODAY      <- format(Sys.time(), "%m%d%Y")

CONTRAST_LABELS <- c( MBP = "MBP vs EtOH")
DEG_COLORS      <- c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey75")

################################################################################
# GENE ANNOTATION
################################################################################

bt <- read.table(GENE_INFO, header = FALSE,
                 col.names = c("gene_id", "gene_name", "gene_type"))
bt$gene_id_clean <- sub("\\..*", "", bt$gene_id)
id_to_name <- setNames(bt$gene_name, bt$gene_id_clean)

################################################################################
# LOAD ALL-GENE RESULTS
# These are the full per-condition result files — all tested genes with stats.
# Significance (adj.P.Val) comes directly from these files; no merge needed.
################################################################################

allgenes <- do.call(rbind, lapply(names(ALLGENE_FILES), function(key) {
  parts <- strsplit(key, "_")[[1]]
  res   <- read.table(ALLGENE_FILES[[key]], header = TRUE, sep = "\t", row.names = 1)
  res$Gene      <- rownames(res)
  res$Contrast  <- parts[1]
  res$Timepoint <- parts[2]
  res
}))

allgenes$gene_id_clean <- sub("\\..*", "", allgenes$Gene)
allgenes$gene_name     <- id_to_name[allgenes$gene_id_clean]
allgenes$gene_name     <- ifelse(is.na(allgenes$gene_name), allgenes$Gene, allgenes$gene_name)
allgenes$Contrast      <- factor(allgenes$Contrast,  levels = c("MBP"))
allgenes$Timepoint     <- factor(allgenes$Timepoint, levels = c("T6", "T24"))
allgenes$condition     <- paste0(CONTRAST_LABELS[as.character(allgenes$Contrast)],
                                 " | ", allgenes$Timepoint)

# Color class based on FDR from full results (not filtered DEG file)
allgenes$vcls <- factor(
  ifelse(allgenes$adj.P.Val < FDR_THRESH & allgenes$logFC > 0, "Up",
  ifelse(allgenes$adj.P.Val < FDR_THRESH & allgenes$logFC < 0, "Down", "NS")),
  levels = c("Up", "Down", "NS"))

CONDITIONS <- unique(allgenes$condition)
cat("All-gene data loaded:", nrow(allgenes), "gene-condition rows\n")

################################################################################
# LOAD & ANNOTATE DEGs (used only for bar plot)
################################################################################

deg <- read.table(DEG_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deg$gene_id_clean <- sub("\\..*", "", deg$Gene)
deg$gene_name     <- id_to_name[deg$gene_id_clean]
deg$gene_name     <- ifelse(is.na(deg$gene_name), deg$Gene, deg$gene_name)

deg$class <- factor(
  ifelse(deg$adj.P.Val < FDR_THRESH & deg$logFC >  FC_THRESH, "Up",
  ifelse(deg$adj.P.Val < FDR_THRESH & deg$logFC < -FC_THRESH, "Down", "NS")),
  levels = c("Up", "Down", "NS"))

deg$Contrast  <- factor(deg$Contrast,  levels = c("MBP"))
deg$Timepoint <- factor(deg$Timepoint, levels = c("T6", "T24"))
deg$condition <- paste0(CONTRAST_LABELS[as.character(deg$Contrast)], " | ", deg$Timepoint)

cat("DEG counts per condition:\n")
for (cond in CONDITIONS) {
  s <- deg[deg$condition == cond, ]
  cat(sprintf("  %-30s  Up: %4d  Down: %4d\n", cond,
              sum(s$class == "Up"), sum(s$class == "Down")))
}

################################################################################
# 1. BAR PLOT
################################################################################

bar_data <- deg %>%
  group_by(condition, Contrast, Timepoint) %>%
  summarise(Up   =  sum(class == "Up"),
            Down = -sum(class == "Down"),
            .groups = "drop") %>%
  pivot_longer(c(Up, Down), names_to = "direction", values_to = "count") %>%
  mutate(direction = factor(direction, levels = c("Up", "Down")))

p_bar <- ggplot(bar_data, aes(x = condition, y = count, fill = direction)) +
  geom_col(width = 0.6, position = "identity") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = abs(count),
                vjust = ifelse(count >= 0, -0.3, 1.2)),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(Up = "#D7191C", Down = "#2C7BB6"),
                    labels = c("Up-regulated", "Down-regulated")) +
  scale_y_continuous(labels = abs) +
  labs(
    title = sprintf("DEG Counts (FDR < %d%%, |log2FC| > %.2f)",
                    FDR_THRESH * 100, FC_THRESH),
    x = NULL, y = "Number of DEGs", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(OUT_DIR, paste0(TODAY, "_barplot_DEGs.png")),
       p_bar, width = 7, height = 5, dpi = 300)
cat("Bar plot saved.\n")

################################################################################
# 2 & 3. VOLCANO + MA — one page per condition
# Volcano plots ALL tested genes; color = FDR significance from full results
################################################################################

for (cond in CONDITIONS) {
  bg <- allgenes[allgenes$condition == cond, ]

  n_up   <- sum(bg$vcls == "Up")
  n_dn   <- sum(bg$vcls == "Down")
  n_up_fc <- sum(bg$adj.P.Val < FDR_THRESH & bg$logFC >  FC_LINE)
  n_dn_fc <- sum(bg$adj.P.Val < FDR_THRESH & bg$logFC < -FC_LINE)
  xrange  <- range(bg$logFC, na.rm = TRUE)
  ymax    <- max(-log10(bg$P.Value), na.rm = TRUE)

  # Label top DEGs by FDR then effect size
  label_df <- bg %>%
    filter(vcls != "NS") %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    slice_head(n = N_LABEL)

  # --- Volcano ---
  pv <- ggplot(bg, aes(x = logFC, y = -log10(P.Value), color = vcls)) +
    geom_point(size = 0.9, alpha = 0.55) +
    geom_text_repel(data = label_df, aes(label = gene_name),
                    size = 2.5, max.overlaps = 25, segment.size = 0.3,
                    show.legend = FALSE) +
    geom_vline(xintercept = c(-FC_LINE, FC_LINE),
               linetype = "solid", linewidth = 0.5, color = "#E69F00") +
    scale_color_manual(values = DEG_COLORS,
                       guide = guide_legend(override.aes = list(size = 2.5))) +
    annotate("text", x = xrange[2], y = ymax,
             label = sprintf("All up: %d\n|FC|>%.1f: %d", n_up, FC_LINE, n_up_fc),
             hjust = 1, vjust = 1, size = 3, color = "#D7191C", fontface = "bold") +
    annotate("text", x = xrange[1], y = ymax,
             label = sprintf("All dn: %d\n|FC|>%.1f: %d", n_dn, FC_LINE, n_dn_fc),
             hjust = 0, vjust = 1, size = 3, color = "#2C7BB6", fontface = "bold") +
    labs(title = paste("Volcano —", cond),
         x = expression(log[2]~FC),
         y = expression(-log[10]~italic(P)),
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  # --- MA ---
  pm <- ggplot(bg, aes(x = AveExpr, y = logFC, color = vcls)) +
    geom_point(size = 0.9, alpha = 0.55) +
    geom_text_repel(data = label_df, aes(label = gene_name),
                    size = 2.5, max.overlaps = 25, segment.size = 0.3,
                    show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_hline(yintercept = c(-FC_THRESH, FC_THRESH),
               linetype = "dashed", linewidth = 0.4, color = "grey40") +
    scale_color_manual(values = DEG_COLORS,
                       guide = guide_legend(override.aes = list(size = 2.5))) +
    labs(title = paste("MA —", cond),
         x = "Average Expression",
         y = expression(log[2]~FC),
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  fname <- gsub("[^A-Za-z0-9_]", "_", cond)
  ggsave(file.path(OUT_DIR, paste0(TODAY, "_volcano_MA_", fname, ".png")),
         pv + pm, width = 13, height = 6, dpi = 300)
}

cat("Volcano + MA plots saved.\n")
cat("Done. Figures in:", OUT_DIR, "\n")
