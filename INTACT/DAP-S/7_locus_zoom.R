#!/usr/bin/env Rscript
# 7_locus_zoom.R
# LocusZoom-style fine-mapping plots for top genes per condition.
#
# Three-panel layout per gene:
#   1. PIP panel      — per-variant PIP, colored by signal cluster
#   2. Chromatin panel — open chromatin regions in the locus
#   3. Gene track      — nearby genes from Gencode v31, strand-aware
#
# Selects top gene per condition by max PIP, outputs one PNG per condition.

library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

###############################################################################
# CONFIG
###############################################################################

DAPSDIR    <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/daps_output"
OUTDIR     <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/daps_output/figures/locus_zoom"
CHROMATIN  <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/TORUS/open_chromatin_hg38_clean.bed"
GENCODE    <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz"
EQTL_DIR   <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_input/split"
WINDOW_KB  <- 100   # cis-window in kb (for display)
PAD_KB     <- 10    # extra padding beyond outermost tested SNP

CONDITIONS <- c("BPA_T6", "BPA_T24", "EtOH_T6", "EtOH_T24",
                "H2O_T6", "H2O_T24", "MBP_T6",  "MBP_T24")

# Signal cluster color palette (up to 8 clusters)
SC_COLS <- c("1"="#D62728", "2"="#123B7A", "3"="#2D6A2D", "4"="#8B4DAB",
             "5"="#D4A017", "6"="#E8649A", "7"="#56B4E9", "8"="#80C840",
             "none"="grey75")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Load Gencode gene annotations (do once)
###############################################################################

cat("Loading Gencode v31 gene annotations...\n")
anno <- fread(
  cmd = paste("zcat", GENCODE, "| grep -v '^#' | awk '$3==\"gene\"'"),
  header = FALSE, sep = "\t",
  col.names = c("chr","source","feature","start","end","score","strand","frame","attributes")
)
anno[, gene_id   := sub("ID=([^;]+).*",     "\\1", attributes)]
anno[, gene_name := sub(".*gene_name=([^;]+).*", "\\1", attributes)]
anno[, gene_base := sub("\\..*", "", gene_id)]
anno[, chr_clean := sub("chr", "", chr)]
anno <- anno[!grepl("PAR_Y", gene_id)]
cat("  Loaded", nrow(anno), "gene records\n\n")

###############################################################################
# Helper: build one locus zoom plot for a given gene + condition
###############################################################################

plot_locus <- function(cond, gene_id, gene_name, pip_dat, eqtl_dat) {

  # --- Locus boundaries ---
  pip_sub  <- pip_dat[gene == gene_id]
  eqtl_sub <- eqtl_dat[gene == gene_id]

  # Parse positions
  pip_sub[, c("chr","pos","ref","alt") := tstrsplit(SNP_ID, ":", fixed=TRUE)]
  pip_sub[, pos := as.integer(pos)]
  chr_num <- pip_sub$chr[1]
  chr_str <- paste0("chr", chr_num)

  xmin <- min(pip_sub$pos) - PAD_KB * 1000
  xmax <- max(pip_sub$pos) + PAD_KB * 1000
  xmid <- (xmin + xmax) / 2

  # Signal cluster as factor for color
  pip_sub[, sc_col := ifelse(is.na(SC_ID), "none", as.character(SC_ID))]
  pip_sub[, sc_col := factor(sc_col, levels = c(as.character(1:8), "none"))]
  pip_sub[, pt_size := ifelse(is.na(SC_ID), 0.8, 2.2)]

  # --------------------------------------------------------------------------
  # Panel 1: PIP
  # --------------------------------------------------------------------------
  p_pip <- ggplot(pip_sub, aes(x = pos / 1e6, y = PIP, color = sc_col,
                                size = pt_size)) +
    geom_segment(aes(xend = pos / 1e6, yend = 0),
                 linewidth = 0.3, alpha = 0.5) +
    geom_point() +
    scale_color_manual(values = SC_COLS, drop = TRUE,
                       name = "Signal\ncluster",
                       na.value = "grey75") +
    scale_size_identity() +
    scale_x_continuous(limits = c(xmin, xmax) / 1e6,
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous("PIP", limits = c(0, 1.05),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    ggtitle(bquote(italic(.(gene_name)) ~ "—" ~ .(cond)),
            subtitle = sprintf("chr%s:%s–%s Mb  |  max PIP = %.3f",
                               chr_num,
                               format(round(xmin/1e6, 3), nsmall=3),
                               format(round(xmax/1e6, 3), nsmall=3),
                               max(pip_sub$PIP))) +
    theme_bw(base_size = 10) +
    theme(axis.title.x    = element_blank(),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title      = element_text(size = 11, face = "bold"),
          plot.subtitle   = element_text(size = 8, color = "grey40"),
          legend.key.size = unit(0.6, "lines"),
          legend.text     = element_text(size = 8),
          legend.title    = element_text(size = 8))

  # --------------------------------------------------------------------------
  # Panel 2: Open chromatin track
  # --------------------------------------------------------------------------
  chrom <- fread(CHROMATIN, header = FALSE, sep = "\t",
                 select = 1:3, col.names = c("chr","start","end"))
  chrom_sub <- chrom[chr == chr_str & start >= xmin & end <= xmax]

  if (nrow(chrom_sub) > 0) {
    p_chrom <- ggplot(chrom_sub) +
      geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6,
                    ymin = 0, ymax = 1),
                fill = "#0F6E56", alpha = 0.7, color = NA) +
      scale_x_continuous(limits = c(xmin, xmax) / 1e6,
                         expand = expansion(mult = 0.01)) +
      scale_y_continuous("ATAC", breaks = NULL) +
      theme_bw(base_size = 10) +
      theme(axis.title.x   = element_blank(),
            axis.text.x    = element_blank(),
            axis.ticks.x   = element_blank(),
            panel.grid      = element_blank(),
            axis.title.y    = element_text(size = 8, angle = 90))
  } else {
    p_chrom <- ggplot() +
      annotate("text", x = xmid / 1e6, y = 0.5,
               label = "No open chromatin in window", size = 3,
               color = "grey50") +
      scale_x_continuous(limits = c(xmin, xmax) / 1e6) +
      scale_y_continuous("ATAC", breaks = NULL) +
      theme_bw(base_size = 10) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), panel.grid = element_blank())
  }

  # --------------------------------------------------------------------------
  # Panel 3: Gene track
  # --------------------------------------------------------------------------
  genes_in_window <- anno[chr_clean == chr_num &
                            start <= xmax & end >= xmin]

  if (nrow(genes_in_window) > 0) {
    # Assign y-levels to avoid overlap (simple stacking)
    genes_in_window <- genes_in_window[order(start)]
    genes_in_window[, y := 0]
    for (j in seq_len(nrow(genes_in_window))) {
      if (j == 1) next
      prev <- genes_in_window[seq_len(j-1)]
      overlap <- prev[end >= genes_in_window$start[j] &
                        y == genes_in_window$y[j]]
      if (nrow(overlap) > 0) genes_in_window[j, y := 1]
    }

    # Focal gene highlighted
    focal_id <- gene_id  # capture in local scope for data.table
    genes_in_window[, is_focal := gene_base == focal_id]
    genes_in_window[, label := gene_name]  # explicit copy avoids scoping issues

    p_gene <- ggplot(genes_in_window) +
      # Gene body
      geom_rect(aes(xmin = pmax(start, xmin) / 1e6,
                    xmax = pmin(end, xmax) / 1e6,
                    ymin = y - 0.25, ymax = y + 0.25,
                    fill = is_focal),
                color = NA, alpha = 0.85) +
      # Strand arrows
      geom_text(aes(x = (pmax(start, xmin) + pmin(end, xmax)) / 2 / 1e6,
                    y = y,
                    label = ifelse(strand == "+", ">", "<")),
                size = 2.5, color = "white") +
      # Gene name label — only for genes wider than 2% of window
      geom_text(
        data = ~ .x[pmin(end, xmax) - pmax(start, xmin) > (xmax - xmin) * 0.02],
        aes(x = (pmax(start, xmin) + pmin(end, xmax)) / 2 / 1e6,
            y = y + 0.45,
            label = label,
            fontface = ifelse(is_focal, "bold.italic", "italic")),
        size = 2.5, color = "grey20", vjust = 0) +
      scale_fill_manual(values = c("TRUE" = "#E24B4A", "FALSE" = "#B4B2A9"),
                        guide = "none") +
      scale_x_continuous(
        sprintf("Genomic position — chr%s (Mb)", chr_num),
        limits = c(xmin, xmax) / 1e6,
        expand = expansion(mult = 0.01)
      ) +
      scale_y_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
      theme_bw(base_size = 10) +
      theme(axis.title.y   = element_blank(),
            panel.grid     = element_blank(),
            axis.title.x   = element_text(size = 9))
  } else {
    p_gene <- ggplot() +
      scale_x_continuous(sprintf("Genomic position — chr%s (Mb)", chr_num),
                         limits = c(xmin, xmax) / 1e6) +
      annotate("text", x = xmid / 1e6, y = 0.5,
               label = "No genes in window", size = 3, color = "grey50") +
      theme_bw(base_size = 10) +
      theme(panel.grid = element_blank(), axis.title.y = element_blank())
  }

  # --------------------------------------------------------------------------
  # Combine panels
  # --------------------------------------------------------------------------
  plot_grid(p_pip, p_chrom, p_gene,
            ncol = 1, align = "v", axis = "lr",
            rel_heights = c(3, 0.9, 1.2))
}

###############################################################################
# Main: select top gene per condition and plot
###############################################################################

cat("Selecting top gene per condition and plotting...\n\n")

for (cond in CONDITIONS) {

  cat("Condition:", cond, "")

  # Load PIP data for this condition
  pip_file <- file.path(DAPSDIR, paste0(cond, "_daps_pips.txt.gz"))
  pip <- fread(pip_file, header = TRUE)
  setnames(pip, c("gene","SC_ID","SNP_ID","PIP"))

  # Load eQTL data for z-scores (not used in plot but used for gene selection)
  eqtl <- fread(file.path(EQTL_DIR, paste0(cond, ".eQTL.txt.gz")),
                header = TRUE, fill = TRUE)
  setnames(eqtl, c("SNP","gene","beta","t.stat","p.value"))

  # Select gene with highest single-variant PIP
  top <- pip[which.max(PIP)]
  top_gene <- top$gene
  top_name <- anno[gene_base == top_gene, gene_name][1]
  if (is.na(top_name)) top_name <- top_gene

  cat("->", top_name, "| max PIP:", round(top$PIP, 3), "\n")

  # Build plot
  p <- tryCatch(
    plot_locus(cond, top_gene, top_name, pip, eqtl),
    error = function(e) {
      cat("  Plot error:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(p)) {
    figfn <- file.path(OUTDIR, paste0("locus_", cond, "_", top_name, ".png"))
    ggsave(figfn, p, width = 7, height = 6, dpi = 150)
    cat("  Written:", basename(figfn), "\n")
  }
}

cat("\nDone.\n")
