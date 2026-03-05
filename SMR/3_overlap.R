#!/usr/bin/env Rscript
# Analyze overlap of SMR with DEGs and reGenes
# Ouputs:
#   1. SNP source bar charts (whether top SNP from treatment or control)
#   2. Effect size comparisons 
#   3. Upset plot of overlap across conditions
#   4. Barplot of overlaps 
#   5. P-value distribution faceted by category

library(tidyverse)
library(data.table)
library(annotables)
library(openxlsx)
library(ggrepel)
library(ComplexHeatmap)
library(UpSetR)
rm(list=ls())

base_dir   <- "/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826"
trait      <- "CARDIoGRAM_C4D_CAD"
fdr_thresh <- 0.1

pairs <- list(
  BPA_T6  = c(treatment="BPA_T6",  control="H2O_T6",  contrast="BPA", tp="T6"),
  BPA_T24 = c(treatment="BPA_T24", control="H2O_T24", contrast="BPA", tp="T24"),
  MBP_T6  = c(treatment="MBP_T6",  control="EtOH_T6", contrast="MBP", tp="T6"),
  MBP_T24 = c(treatment="MBP_T24", control="EtOH_T24",contrast="MBP", tp="T24")
)

# eGene removed from all categories
cat_colors <- c(
  "DEG + reGene" = "#6a3d9a",
  "reGene"       = "#ff7f00",
  "DEG"          = "#1f78b4",
  "TWAS only"    = "#7d7a7a"
)

outdir <- "./3_overlap_analysis/"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

grch38_unq <- grch38 %>%
  filter(chr %in% as.character(1:22), grepl("protein", biotype)) %>%
  distinct(ensgene, .keep_all=TRUE) %>%
  dplyr::select(gene=ensgene, symbol)

###############################################################################
### Load supporting data
###############################################################################

cat("Loading mashr SNP classification...\n")
snp_class <- fread(paste0("./1_SMR_output/", trait, "/", trait, "_snp_classification.txt.gz"),
                   header=TRUE, data.table=FALSE)

cat("Loading DEGs...\n")
deg_all <- fread(file.path(base_dir, "DEG/02222026_AllDEGs_FDR10.txt"), data.table=FALSE)

regenes <- map(names(pairs), function(pair_name) {
  snp_class %>%
    filter(pair == pair_name, gene_class == "response_eGene") %>%
    pull(phenotype_id) %>% unique()
}) %>% setNames(names(pairs))

degs <- map(names(pairs), function(pair_name) {
  p <- pairs[[pair_name]]
  deg_all %>%
    filter(Contrast == p["contrast"], Timepoint == p["tp"],
           adj.P.Val < fdr_thresh) %>%
    pull(Gene) %>% sub("\\..*", "", .) %>% unique()
}) %>% setNames(names(pairs))

annotate_res <- function(res, pair_name) {
  res %>%
    mutate(
      symbol    = ifelse(is.na(symbol) | symbol == "", phenotype_id, symbol),
      is_DEG    = phenotype_id %in% degs[[pair_name]],
      is_reGene = phenotype_id %in% regenes[[pair_name]],
      category  = case_when(
        is_DEG & is_reGene ~ "DEG + reGene",
        is_reGene          ~ "reGene",
        is_DEG             ~ "DEG",
        TRUE               ~ "TWAS only"
      ),
      category = factor(category, levels=names(cat_colors))
    ) %>%
    rename(snp_source = selected_from)
}

load_twas <- function(pair_name, sig_only=TRUE) {
  res <- fread(paste0("./2_summary_twas/", pair_name, "/", trait, "_", pair_name, "_twas.txt"),
               header=TRUE, data.table=FALSE, fill=TRUE)
  if (sig_only) res <- filter(res, FDR < fdr_thresh)
  annotate_res(res, pair_name)
}

###############################################################################
### Part 1: SNP source bar chart + effect size scatter per pair
###############################################################################

for (pair_name in names(pairs)) {
  trt <- pairs[[pair_name]]["treatment"]
  ctl <- pairs[[pair_name]]["control"]
  res <- load_twas(pair_name)

  # 1a: SNP source bar chart
  source_summ <- res %>%
    count(snp_source) %>%
    mutate(pct   = round(100 * n / sum(n), 1),
           label = paste0(n, "\n(", pct, "%)"))

  p_source <- ggplot(source_summ, aes(x=snp_source, y=n, fill=snp_source)) +
    geom_bar(stat="identity", width=0.6) +
    geom_text(aes(label=label), vjust=-0.3, size=4) +
    scale_fill_manual(values=c("treatment"="#e31a1c", "control"="#1f78b4"),
                      labels=c("treatment"=trt, "control"=ctl)) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
    ylab("Number of TWAS genes") +
    ggtitle(paste(pair_name, "— SNP nominated by")) +
    theme_bw(base_size=13) +
    theme(legend.title       = element_blank(),
          axis.title.x       = element_blank(),
          axis.text.x        = element_text(size=11),
          plot.title         = element_text(hjust=0.5, face="bold"),
          panel.grid.major.x = element_blank())

  ggsave(paste0(outdir, "Figure1a_", pair_name, "_snp_source.png"),
         p_source, width=5, height=4.5, dpi=150)

  # 1b: Effect size scatter
  lim <- max(abs(c(res$pm_trt, res$pm_ctl)), na.rm=TRUE) * 1.15

  cat_counts   <- res %>% count(category) %>%
    mutate(legend_label=paste0(category, "  (n=", n, ")"))
  color_labels <- setNames(cat_counts$legend_label, as.character(cat_counts$category))

  label_dat <- res %>%
    mutate(
      n_cats      = is_DEG + is_reGene,
      dist_origin = sqrt(pm_trt^2 + pm_ctl^2),
      dist_diag   = abs(pm_trt - pm_ctl) / sqrt(2),
      score       = scale(dist_origin)[,1] * 0.5 +
                    scale(dist_diag)[,1]   * 0.35 +
                    scale(n_cats)[,1]      * 0.15
    ) %>%
    arrange(desc(score)) %>%
    slice_head(n=10)

  p_scatter <- ggplot(res %>% arrange(is_DEG + is_reGene),
                      aes(x=pm_ctl, y=pm_trt, color=category, shape=snp_source)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey70", linewidth=0.5) +
    geom_hline(yintercept=0, color="grey85", linewidth=0.4) +
    geom_vline(xintercept=0, color="grey85", linewidth=0.4) +
    geom_point(data=~filter(., category=="TWAS only"), alpha=0.4, size=1.8) +
    geom_point(data=~filter(., category!="TWAS only"), alpha=0.9, size=2.5) +
    geom_text_repel(data=label_dat, aes(label=symbol),
                    size=3, fontface="italic",
                    box.padding=0.6, point.padding=0.3,
                    segment.color="grey50", segment.size=0.3,
                    max.overlaps=Inf, show.legend=FALSE) +
    scale_color_manual(values=cat_colors, labels=color_labels, drop=TRUE) +
    scale_shape_manual(values=c("treatment"=17, "control"=16),
                       labels=c("treatment"=paste("SNP from", trt),
                                "control"=paste("SNP from", ctl))) +
    coord_fixed(xlim=c(-lim, lim), ylim=c(-lim, lim)) +
    xlab(paste("Posterior mean effect —", ctl)) +
    ylab(paste("Posterior mean effect —", trt)) +
    ggtitle(pair_name) +
    guides(color=guide_legend(order=1, override.aes=list(size=3)),
           shape=guide_legend(order=2, override.aes=list(size=3))) +
    theme_bw(base_size=13) +
    theme(plot.title       = element_text(hjust=0.5, size=14, face="bold"),
          legend.title     = element_blank(),
          legend.key.size  = unit(0.5, "cm"),
          legend.text      = element_text(size=9),
          panel.grid.major = element_line(color="grey93"),
          panel.grid.minor = element_blank(),
          axis.title       = element_text(size=11))

  ggsave(paste0(outdir, "Figure1b_", pair_name, "_scatter.png"),
         p_scatter, width=8, height=6.5, dpi=150)

  # 1c: Paired dot/line — top 30 by FDR
  plot_df <- res %>%
    slice_min(FDR, n=30, with_ties=FALSE) %>%
    dplyr::select(phenotype_id, symbol, pm_trt, pm_ctl) %>%
    pivot_longer(cols=c(pm_trt, pm_ctl), names_to="condition", values_to="effect") %>%
    mutate(condition=recode(condition, pm_trt=trt, pm_ctl=ctl))

  p_paired <- ggplot(plot_df, aes(x=condition, y=effect, group=phenotype_id)) +
    geom_line(alpha=0.4, color="grey50") +
    geom_point(aes(color=condition), size=2) +
    geom_hline(yintercept=0, linetype="dashed", color="grey60") +
    ylab("Posterior mean eQTL effect") +
    ggtitle(paste(pair_name, "— top 30 genes by FDR")) +
    theme_bw(base_size=13) +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          plot.title      = element_text(hjust=0.5, size=10))

  ggsave(paste0(outdir, "Figure1c_", pair_name, "_paired.png"),
         p_paired, width=4, height=5, dpi=150)

  cat("Plots saved:", pair_name, "\n")
}

###############################################################################
### Part 2: Cross-pair overlap — UpSet + binary heatmap
###############################################################################

twas <- map(names(pairs), function(p) load_twas(p) %>% pull(phenotype_id)) %>%
        setNames(names(pairs))

all_genes <- unique(unlist(twas))
bin_mat   <- sapply(names(pairs), function(p) as.integer(all_genes %in% twas[[p]]))
rownames(bin_mat) <- all_genes

# 2a: UpSet
png(paste0(outdir, "Figure2a_upset.png"), width=800, height=500, res=120)
print(upset(as.data.frame(bin_mat), sets=names(pairs), order.by="freq",
            mainbar.y.label="Gene intersections", sets.x.label="TWAS genes per pair",
            text.scale=1.2))
dev.off()
cat("UpSet plot saved\n")

# 2b: Binary heatmap — genes in ≥2 pairs
shared_genes <- all_genes[rowSums(bin_mat) >= 2]
if (length(shared_genes) > 0) {
  mat_shared <- bin_mat[shared_genes, ]
  sym_map    <- setNames(grch38_unq$symbol, grch38_unq$gene)
  row_labels <- ifelse(is.na(sym_map[shared_genes]) | sym_map[shared_genes] == "",
                       shared_genes, sym_map[shared_genes])

  ht <- Heatmap(mat_shared,
                col=c("0"="white", "1"="#2171b5"),
                name="TWAS sig",
                cluster_rows=TRUE, cluster_columns=FALSE,
                show_row_names=nrow(mat_shared) <= 80,
                row_labels=row_labels,
                row_names_gp=gpar(fontsize=6),
                column_names_gp=gpar(fontsize=10),
                heatmap_legend_param=list(
                  at=c(0,1), labels=c("No","Yes"),
                  title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=8)))

  png(paste0(outdir, "Figure2b_heatmap_overlap.png"),
      width=500, height=max(400, length(shared_genes) * 8), res=120)
  draw(ht)
  dev.off()
  cat("Binary heatmap saved:", length(shared_genes), "genes in ≥2 pairs\n")
}

###############################################################################
### Part 3: TWAS overlap with DEG / reGene — bar chart + summary tables
###############################################################################

overlap_df <- map_dfr(names(pairs), function(pair_name) {
  load_twas(pair_name) %>%
    dplyr::select(pair, gene=phenotype_id, is_DEG, is_reGene, category)
})

write.xlsx(overlap_df,   paste0(outdir, "1_twas_overlap_categories.xlsx"), overwrite=TRUE)

summ_overlap <- overlap_df %>%
  group_by(pair) %>%
  summarize(n_twas   = n(),
            n_DEG    = sum(is_DEG),
            n_reGene = sum(is_reGene),
            n_both   = sum(is_DEG & is_reGene),
            .groups  = "drop")

write.xlsx(summ_overlap, paste0(outdir, "2_twas_overlap_summary.xlsx"), overwrite=TRUE)
print(summ_overlap)

plot_summ <- overlap_df %>%
  count(pair, category) %>%
  mutate(category=factor(category, levels=names(cat_colors)))

p_bar <- ggplot(plot_summ, aes(x=pair, y=n, fill=category)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=n), position=position_stack(vjust=0.5), size=2.8) +
  scale_fill_manual(values=cat_colors) +
  ylab("Number of TWAS genes (FDR<10%)") +
  ggtitle("TWAS gene overlap with DEG and reGene") +
  theme_bw(base_size=13) +
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_text(angle=45, hjust=1),
        legend.title    = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        plot.title      = element_text(hjust=0.5, size=11))

ggsave(paste0(outdir, "Figure3_overlap_barplot.png"),
       p_bar, width=7, height=5, dpi=150)

###############################################################################
### Part 4: P-value distribution faceted by category (all genes, no FDR filter)
# Vertical lines: dashed grey = nominal p=0.05, dotted red = FDR 10% threshold
###############################################################################

pval_df <- map_dfr(names(pairs), function(pair_name) {
  load_twas(pair_name, sig_only=FALSE) %>%
    dplyr::select(pair, phenotype_id, pval_gwas, FDR, category)
}) %>% filter(!is.na(pval_gwas))

# Per facet: p-value at FDR 10% boundary (max p among FDR-significant genes)
fdr_lines <- pval_df %>%
  group_by(pair, category) %>%
  summarize(
    fdr_pval = {
      sig <- pval_gwas[!is.na(FDR) & FDR < fdr_thresh]
      if (length(sig) > 0) max(sig) else NA_real_
    },
    .groups = "drop"
  ) %>%
  filter(!is.na(fdr_pval))

p_pvaldist <- ggplot(pval_df, aes(x=pval_gwas)) +
  geom_histogram(aes(fill=category), bins=40, color=NA, alpha=0.85) +
  geom_vline(xintercept=0.05, linetype="dashed", color="grey40", linewidth=0.6) +
  geom_vline(data=fdr_lines, aes(xintercept=fdr_pval),
             linetype="dotted", color="#e31a1c", linewidth=0.6) +
  scale_fill_manual(values=cat_colors) +
  scale_x_continuous("TWAS p-value", expand=c(0, 0)) +
  scale_y_continuous("Number of genes", expand=expansion(mult=c(0, 0.15))) +
  facet_grid(category ~ pair, scales="free_y") +
  theme_bw(base_size=11) +
  theme(legend.position = "none",
        strip.text      = element_text(size=9),
        axis.text.x     = element_text(size=8, angle=45, hjust=1),
        panel.grid      = element_blank(),
        plot.title      = element_text(hjust=0.5, face="bold")) +
  ggtitle("TWAS p-value distribution by category\n(grey dashed = p<0.05 | red dotted = FDR<10%)")

ggsave(paste0(outdir, "Figure4_pval_distribution.png"),
       p_pvaldist,
       width  = 3 * length(pairs),
       height = 2.5 * length(cat_colors),
       dpi    = 150)

cat("P-value distribution plot saved\n")
cat("All done. Outputs in", outdir, "\n")
