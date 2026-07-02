#!/usr/bin/env Rscript
################################################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(data.table)
  library(annotables)
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript GO_ORR.R <CONTRAST> <TIMEPOINT>")

deg_contrast <- args[1]
timepoint    <- args[2]


################################################################################
# CONFIG
################################################################################

BASE_DIR  <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates"
DEG_DIR   <- file.path(BASE_DIR, "DEG")
OUT_DIR     <- file.path(DEG_DIR, "GO_enrichment")
FIGURES_DIR <- file.path(DEG_DIR, "figures")
dir.create(OUT_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE)

TITLE     <- paste0(deg_contrast, "_T", timepoint)
FDR       <- 0.10
TOP_N     <- 15   # max terms per facet
TODAY     <- format(Sys.time(), "%m%d%Y")

cat("\n=== GO/KEGG ORA:", TITLE, "===\n")

################################################################################
# ANNOTATIONS
################################################################################

anno <- grch38 %>%
  dplyr::select(ensgene, entrez) %>%
  distinct() %>%
  filter(!is.na(entrez))

entrez_vec <- function(ensgenes) as.character(filter(anno, ensgene %in% ensgenes)$entrez)

################################################################################
# LOAD DATA
################################################################################

# DEGs — full results file as background; split sig genes by direction
deg_file <- tail(sort(list.files(DEG_DIR,
  pattern = paste0(".*_", deg_contrast, "_T", timepoint, ".txt"),
  full.names = TRUE)), 1)
cat("DEG file:", basename(deg_file), "\n")

degs <- fread(deg_file)
if (colnames(degs)[1] %in% c("V1", "", "Gene", "gene", "gene_id")) {
  setnames(degs, 1, "ensgene")
}
degs <- degs %>%
  filter(!is.na(adj.P.Val)) %>%
  inner_join(anno, by = "ensgene")

deg_background <- as.character(degs$entrez)
deg_up   <- degs %>% filter(adj.P.Val < FDR, logFC >  0) %>% pull(entrez) %>% as.character()
deg_down <- degs %>% filter(adj.P.Val < FDR, logFC <  0) %>% pull(entrez) %>% as.character()
deg_all  <- degs %>% filter(adj.P.Val < FDR)              %>% pull(entrez) %>% as.character()
cat(sprintf("DEGs — Up: %d  Down: %d  All: %d  Background: %d\n",
            length(deg_up), length(deg_down), length(deg_all), length(deg_background)))

 
################################################################################
# ENRICHMENT HELPERS
################################################################################
 
# Run GO BP ORA; returns tidy data.frame or NULL
run_go <- function(genes, bg) {
  if (length(genes) < 5) return(NULL)
  res <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "BP",
                  universe = bg, readable = TRUE)
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  res <- clusterProfiler::simplify(res, cutoff = 0.6, by = "p.adjust")
  as_tibble(res@result) %>%
    filter(qvalue < FDR) %>%
    mutate(
      GeneRatio  = as.numeric(sub("/.*", "", GeneRatio)) /
                   as.numeric(sub(".*/", "", GeneRatio)),
      BgRatio    = as.numeric(sub("/.*", "", BgRatio)) /
                   as.numeric(sub(".*/", "", BgRatio)),
      FoldEnrich = GeneRatio / BgRatio
    )
}
 
# Run KEGG ORA; returns tidy data.frame or NULL
run_kegg <- function(genes, bg) {
  if (length(genes) < 5) return(NULL)
  res <- enrichKEGG(gene = genes, organism = "hsa",
                    minGSSize = 5, universe = bg)
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  as_tibble(res@result) %>%
    filter(qvalue < FDR) %>%
    mutate(
      GeneRatio  = as.numeric(sub("/.*", "", GeneRatio)) /
                   as.numeric(sub(".*/", "", GeneRatio)),
      BgRatio    = as.numeric(sub("/.*", "", BgRatio)) /
                   as.numeric(sub(".*/", "", BgRatio)),
      FoldEnrich = GeneRatio / BgRatio
    )
}
 
################################################################################
# COLOR SCALES
# Up-regulated / All DEGs: light pink -> dark red
# Down-regulated:          light blue -> dark blue
################################################################################
 
color_scale_up   <- scale_color_gradient(low = "#fee0d2", high = "#a50f15",
                                         name = expression(-log[10](q)))
color_scale_down <- scale_color_gradient(low = "#deebf7", high = "#08306b",
                                         name = expression(-log[10](q)))
 
################################################################################
# PLOT HELPER
################################################################################
 
dot_plot <- function(df, title, facet_col = NULL) {
  df <- df %>%
    group_by(across(any_of(facet_col))) %>%
    slice_min(order_by = p.adjust, n = TOP_N) %>%
    ungroup() %>%
    mutate(
      Description = str_wrap(Description, width = 45),
      neg_log10_q = -log10(qvalue + 1e-10)
    )
 
  # Order terms so most significant (lowest qvalue) appears at top of y-axis
  if (!is.null(facet_col)) {
    df <- df %>%
      mutate(Description = factor(Description,
               levels = df %>%
                 arrange(get(facet_col), qvalue) %>%
                 pull(Description) %>% unique()))
  } else {
    df <- df %>% mutate(Description = fct_reorder(Description, desc(qvalue)))
  }
 
  # Separate plots per direction assembled with patchwork
  if (!is.null(facet_col) && facet_col == "Direction") {
    dirs   <- levels(df$Direction)
    plots  <- lapply(dirs, function(d) {
      sub_df <- filter(df, Direction == d)
      if (nrow(sub_df) == 0) return(NULL)
      clr <- if (d == "Down-regulated") color_scale_down else color_scale_up
      ggplot(sub_df, aes(x = FoldEnrich, y = Description,
                         size = Count, color = neg_log10_q)) +
        geom_point() +
        clr +
        scale_size_continuous(name = "Gene count", range = c(2, 8)) +
        labs(title = d, x = "Fold enrichment", y = NULL) +
        theme_bw(base_size = 11) +
        theme(
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border       = element_rect(color = "black"),
          axis.text          = element_text(color = "black"),
          axis.title         = element_text(color = "black"),
          axis.ticks         = element_line(color = "black"),
          strip.background   = element_rect(fill = "grey92", color = "black"),
          strip.text         = element_text(face = "bold", color = "black"),
          plot.title         = element_text(face = "bold", size = 11,
                                            color = "black"),
          legend.text        = element_text(color = "black"),
          legend.title       = element_text(color = "black")
        )
    })
    plots <- Filter(Negate(is.null), plots)
    return(patchwork::wrap_plots(plots, nrow = 1) +
             patchwork::plot_annotation(title = title,
               theme = theme(plot.title = element_text(face = "bold",
                                                       color = "black"))))
  }
 
  # Non-directional plot (All DEGs or single group)
  p <- ggplot(df, aes(x = FoldEnrich, y = Description,
                      size = Count, color = neg_log10_q)) +
    geom_point() +
    color_scale_up +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    labs(title = title, x = "Fold enrichment", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border       = element_rect(color = "black"),
      axis.text          = element_text(color = "black"),
      axis.title         = element_text(color = "black"),
      axis.ticks         = element_line(color = "black"),
      strip.background   = element_rect(fill = "grey92", color = "black"),
      strip.text         = element_text(face = "bold", color = "black"),
      plot.title         = element_text(face = "bold", size = 12,
                                        color = "black"),
      legend.text        = element_text(color = "black"),
      legend.title       = element_text(color = "black")
    )
  p
}
 
################################################################################
# SAVE HELPERS
################################################################################
 
save_results <- function(df, label, type) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  out <- file.path(OUT_DIR, paste0(TODAY, "_", label, "_", type, ".txt"))
  write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", basename(out), "\n")
}
 
save_plot <- function(p, label, type, w = 10, h = 7) {
  out <- file.path(FIGURES_DIR, paste0(TODAY, "_", label, "_", type, ".png"))
  ggsave(out, p, width = w, height = h, dpi = 300)
  cat("Saved:", basename(out), "\n")
}
 
################################################################################
# DEG ENRICHMENT — split Up/Down/All, direction-aware colors
################################################################################
 
cat("\n--- DEG enrichment ---\n")
 
# GO BP
go_up   <- run_go(deg_up,   deg_background)
go_down <- run_go(deg_down, deg_background)
go_all  <- run_go(deg_all,  deg_background)
 
save_results(go_up,   paste0(TITLE, "_DEGs_Up"),   "GO_BP")
save_results(go_down, paste0(TITLE, "_DEGs_Down"), "GO_BP")
save_results(go_all,  paste0(TITLE, "_DEGs_All"),  "GO_BP")
 
go_deg_combined <- bind_rows(
  mutate(go_up,   Direction = "Up-regulated"),
  mutate(go_down, Direction = "Down-regulated"),
  mutate(go_all,  Direction = "All DEGs")
) %>% mutate(Direction = factor(Direction,
               levels = c("Up-regulated", "Down-regulated", "All DEGs")))
 
if (nrow(go_deg_combined) > 0) {
  p <- dot_plot(go_deg_combined,
                title = paste(TITLE, "— GO BP (ORA, FDR < 10%)"),
                facet_col = "Direction")
  save_plot(p, paste0(TITLE, "_DEGs"), "GO_BP_dotplot", w = 18, h = 7)
} else {
  cat("No significant GO BP terms for DEGs\n")
}
 
# KEGG
kegg_up   <- run_kegg(deg_up,   deg_background)
kegg_down <- run_kegg(deg_down, deg_background)
kegg_all  <- run_kegg(deg_all,  deg_background)
 
save_results(kegg_up,   paste0(TITLE, "_DEGs_Up"),   "KEGG")
save_results(kegg_down, paste0(TITLE, "_DEGs_Down"), "KEGG")
save_results(kegg_all,  paste0(TITLE, "_DEGs_All"),  "KEGG")
 
kegg_deg_combined <- bind_rows(
  mutate(kegg_up,   Direction = "Up-regulated"),
  mutate(kegg_down, Direction = "Down-regulated"),
  mutate(kegg_all,  Direction = "All DEGs")
) %>% mutate(Direction = factor(Direction,
               levels = c("Up-regulated", "Down-regulated", "All DEGs")))
 
if (nrow(kegg_deg_combined) > 0) {
  p <- dot_plot(kegg_deg_combined,
                title = paste(TITLE, "— KEGG (ORA, FDR < 10%)"),
                facet_col = "Direction")
  save_plot(p, paste0(TITLE, "_DEGs"), "KEGG_dotplot", w = 18, h = 7)
} else {
  cat("No significant KEGG pathways for DEGs\n")
}
