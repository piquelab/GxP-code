#!/usr/bin/env Rscript
# TensorQTL Analysis - With Additional FDR Correction on Beta-Approximated P-values
# Generated using Claude 9/9/2025

library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)

# =============================================================================
# CONFIGURATION
# =============================================================================
base_dir <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor"
input_dir <- file.path(base_dir, "tensorqtl_output_cis_SV15_100kb")
output_dir <- file.path(input_dir, "analysis_results_fdr")
figures_dir <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Define conditions
conditions <- c("BPA_100nM_24", "BPA_100nM_6", "MBP_500nM_24", "MBP_500nM_6")

cat("Starting tensorQTL analysis with FDR correction on beta p-values...\n")

# =============================================================================
# 1. LOAD AND PROCESS DATA WITH FDR CORRECTION
# =============================================================================
cat("Loading tensorQTL results...\n")

qtl_data <- list()
qtl_summary <- data.frame()

for (condition in conditions) {
  cat("Processing:", condition, "\n")

  qtl_file <- file.path(input_dir, paste0(condition, "_cis.cis_qtl.txt.gz"))

  if (file.exists(qtl_file)) {
    qtl_results <- fread(qtl_file)
    qtl_results$condition <- condition

    # Check which p-value columns are available
    available_cols <- colnames(qtl_results)
    cat("  Available p-value columns:", paste(grep("pval", available_cols, value = TRUE), collapse = ", "), "\n")
    
    # Use beta-approximated p-values as base for FDR correction, use perm p-val or nominal-pval as back up
    if ("pval_beta" %in% available_cols) {
      qtl_results$pval_base <- qtl_results$pval_beta
      correction_method <- "beta_approximation + FDR"
    } else if ("pval_perm" %in% available_cols) {
      qtl_results$pval_base <- qtl_results$pval_perm
      correction_method <- "permutation + FDR"
    } else {
      cat("WARNING: No corrected p-values found, using nominal p-values + FDR\n")
      qtl_results$pval_base <- qtl_results$pval_nominal
      correction_method <- "nominal + FDR"
    }
    
    # APPLY FDR CORRECTION across all genes in this condition
    cat("Applying FDR correction across", nrow(qtl_results), "genes...\n")
    qtl_results$pval_double_fdr <- p.adjust(qtl_results$pval_base, method = "fdr")
    
    # Also calculate Bonferroni for comparison
    qtl_results$pval_bonferroni <- p.adjust(qtl_results$pval_base, method = "bonferroni")
    
    cat("  Correction method:", correction_method, "\n")
    qtl_data[[condition]] <- qtl_results

    # Calculate summary statistics using FDR values
    n_tests <- nrow(qtl_results)
    n_genes <- length(unique(qtl_results$phenotype_id))
    n_variants <- length(unique(qtl_results$variant_id))

    # Count significant eQTLs using various thresholds
    # Base (single-corrected) thresholds
    base_05 <- qtl_results$pval_base < 0.05 & !is.na(qtl_results$pval_base)
    base_01 <- qtl_results$pval_base < 0.01 & !is.na(qtl_results$pval_base)
    
    #FDR-corrected thresholds
    double_fdr_05 <- qtl_results$pval_double_fdr < 0.05 & !is.na(qtl_results$pval_double_fdr)
    double_fdr_10 <- qtl_results$pval_double_fdr < 0.10 & !is.na(qtl_results$pval_double_fdr)
    double_fdr_20 <- qtl_results$pval_double_fdr < 0.20 & !is.na(qtl_results$pval_double_fdr)
    
    # Bonferroni (for comparison)
    bonf_05 <- qtl_results$pval_bonferroni < 0.05 & !is.na(qtl_results$pval_bonferroni)
    
    # Nominal for comparison
    nom_001 <- qtl_results$pval_nominal < 0.001 & !is.na(qtl_results$pval_nominal)

    # Count eGenes for each threshold
    egenes_base_05 <- length(unique(qtl_results$phenotype_id[base_05]))
    egenes_base_01 <- length(unique(qtl_results$phenotype_id[base_01]))
    egenes_double_fdr_05 <- length(unique(qtl_results$phenotype_id[double_fdr_05]))
    egenes_double_fdr_10 <- length(unique(qtl_results$phenotype_id[double_fdr_10]))
    egenes_double_fdr_20 <- length(unique(qtl_results$phenotype_id[double_fdr_20]))
    egenes_bonf_05 <- length(unique(qtl_results$phenotype_id[bonf_05]))
    egenes_nom_001 <- length(unique(qtl_results$phenotype_id[nom_001]))

    summary_row <- data.frame(
      condition = condition,
      n_tests = n_tests,
      n_genes = n_genes,
      n_variants = n_variants,
      correction_method = correction_method,
      egenes_base_p05 = egenes_base_05,
      egenes_base_p01 = egenes_base_01,
      egenes_double_fdr_05 = egenes_double_fdr_05,
      egenes_double_fdr_10 = egenes_double_fdr_10,
      egenes_double_fdr_20 = egenes_double_fdr_20,
      egenes_bonferroni_05 = egenes_bonf_05,
      egenes_nominal_001 = egenes_nom_001,
      eqtls_base_p05 = sum(base_05),
      eqtls_base_p01 = sum(base_01),
      eqtls_double_fdr_05 = sum(double_fdr_05),
      eqtls_double_fdr_10 = sum(double_fdr_10),
      eqtls_double_fdr_20 = sum(double_fdr_20),
      eqtls_bonferroni_05 = sum(bonf_05),
      eqtls_nominal_001 = sum(nom_001)
    )

    qtl_summary <- rbind(qtl_summary, summary_row)
  } else {
    cat("WARNING: File not found:", qtl_file, "\n")
  }
}

cat("Data loading complete.\n")

# Write FDR summary table
fwrite(qtl_summary, file.path(output_dir, "eQTL_summary_fdr.txt"),
       sep = '\t', quote = FALSE, row.names = FALSE)

sig_threshold <- 0.1
all_sig_genes <- unlist(lapply(qtl_data, function(x) {
  x$phenotype_id[x$pval_double_fdr < sig_threshold] 
}))

# Get unique significant eGenes
unique_sig_egenes <- unique(all_sig_genes)

# Write to file
writeLines(unique_sig_egenes, file.path(input_dir, "unique_significant_egenes.txt"))

# =============================================================================
# 2. QQ PLOTS WITH FDR CORRECTION - may want to use nominal or beta-corrected
# =============================================================================
cat("Creating QQ plots with FDR-corrected p-values...\n")

create_double_fdr_qq_plot <- function(condition_data, condition_name, output_path) {
  # Use FDR-corrected p-values
  pvals <- condition_data$pval_double_fdr[!is.na(condition_data$pval_double_fdr)]
  pvals <- pvals[pvals > 0]
  pvals <- sort(pvals)
  
  n <- length(pvals)
  expected <- (1:n) / n
  pvals[pvals < 1e-20] <- 1e-20
  
  qq_data <- data.frame(
    expected = -log10(expected),
    observed = -log10(pvals)
  )
  
  # Confidence intervals
  ci <- 0.95
  qq_data$clower <- -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n - (1:n) + 1))
  qq_data$cupper <- -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n - (1:n) + 1))
  
  p <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2, fill = "gray") +
    geom_point(alpha = 0.6, size = 0.8, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    theme_classic() +
    labs(title = paste("QQ Plot (FDR-corrected) -", condition_name),
         subtitle = "FDR correction applied to beta-approximated p-values",
         x = expression(Expected -log[10](p)),
         y = expression(Observed -log[10](p))) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"))
  
  ggsave(output_path, p, width = 8, height = 8, dpi = 300)
  return(p)
}

# Create FDR QQ plots
for (condition in conditions) {
  if (condition %in% names(qtl_data)) {  
    output_file <- file.path(figures_dir, paste0("QQ_plot_fdr_", condition, ".png"))
    create_double_fdr_qq_plot(qtl_data[[condition]], condition, output_file)
  }
}

# =============================================================================
# 3. COMPARISON PLOTS: MULTIPLE CORRECTION METHODS
# =============================================================================
cat("Creating multi-method comparison plots...\n")

# Prepare data for comparison
plot_data <- qtl_summary
plot_data$treatment <- sapply(strsplit(plot_data$condition, "_"), function(x) x[1])

# Create comprehensive comparison plot
comparison_long <- plot_data %>%
  select(condition, treatment, egenes_nominal_001, egenes_base_p05, 
         egenes_double_fdr_05, egenes_double_fdr_10, egenes_bonferroni_05) %>%
  pivot_longer(cols = starts_with("egenes_"), 
               names_to = "method", values_to = "egenes") %>%
  mutate(method = case_when(
    method == "egenes_nominal_001" ~ "Nominal p < 0.001",
    method == "egenes_base_p05" ~ "Base corrected p < 0.05", 
    method == "egenes_double_fdr_05" ~ "FDR < 0.05",
    method == "egenes_double_fdr_10" ~ "FDR < 0.10",
    method == "egenes_bonferroni_05" ~ "Bonferroni < 0.05",
    TRUE ~ method
  ))

p_multi_comparison <- ggplot(comparison_long, aes(x = condition, y = egenes, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_viridis_d(name = "Correction\nMethod") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(title = "eGenes Across Multiple Correction Methods",
       subtitle = "Comparison of multiple test corrections",
       x = "Condition", y = "Number of eGenes")

ggsave(file.path(figures_dir, "eGenes_multi_method_comparison.png"), p_multi_comparison,
       width = 14, height = 8, dpi = 300)

# =============================================================================
# 4. P-VALUE DISTRIBUTION COMPARISON
# =============================================================================
cat("Creating p-value distribution comparison...\n")

# Sample one condition for p-value distribution comparison
sample_condition <- names(qtl_data)[1]
sample_data <- qtl_data[[sample_condition]]

p_dist_data <- data.frame(
  base = sample_data$pval_base,
  double_fdr = sample_data$pval_double_fdr,
  bonferroni = sample_data$pval_bonferroni
) %>%
  pivot_longer(everything(), names_to = "method", values_to = "pvalue") %>%
  filter(!is.na(pvalue), pvalue > 0)

p_pval_dist <- ggplot(p_dist_data, aes(x = pvalue, fill = method)) +
  geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
  facet_wrap(~ method, scales = "free_y") +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = paste("P-value Distributions -", sample_condition),
       subtitle = "Comparison of correction methods",
       x = "P-value", y = "Count") +
  scale_x_log10()

ggsave(file.path(figures_dir, paste0("pvalue_distributions_", sample_condition, ".png")), 
       p_pval_dist, width = 12, height = 6, dpi = 300)

# =============================================================================
# 5. FINAL SUMMARY REPORT
# =============================================================================
cat("Generating FDR analysis summary...\n")

if (nrow(qtl_summary) > 0) {
  total_base_05 <- sum(qtl_summary$egenes_base_p05, na.rm = TRUE)
  total_double_fdr_05 <- sum(qtl_summary$egenes_double_fdr_05, na.rm = TRUE)
  total_double_fdr_10 <- sum(qtl_summary$egenes_double_fdr_10, na.rm = TRUE)
  total_bonf_05 <- sum(qtl_summary$egenes_bonferroni_05, na.rm = TRUE)
  total_nominal <- sum(qtl_summary$egenes_nominal_001, na.rm = TRUE)
  total_tests <- sum(qtl_summary$n_tests, na.rm = TRUE)
  
  summary_report <- paste0(
    "FDR CORRECTION TensorQTL Analysis Summary\n",
    paste(rep("=", 60), collapse = ""), "\n",
    "Total conditions analyzed: ", nrow(qtl_summary), "\n",
    "Total statistical tests: ", format(total_tests, big.mark = ","), "\n\n",
    "=== COMPARISON OF CORRECTION METHODS ===\n",
    "Nominal p < 0.001: ", total_nominal, " eGenes\n",
    "Base corrected p < 0.05: ", total_base_05, " eGenes\n",
    "FDR < 0.05: ", total_double_fdr_05, " eGenes\n",
    "FDR < 0.10: ", total_double_fdr_10, " eGenes\n",
    "Bonferroni < 0.05: ", total_bonf_05, " eGenes\n\n",
    "Results saved to: ", output_dir, "\n",
    "Figures saved to: ", figures_dir, "\n"
  )
} else {
  summary_report <- "ERROR: No data could be processed.\n"
}

writeLines(summary_report, file.path(output_dir, "FDR_ANALYSIS_SUMMARY.txt"))
cat(summary_report)

cat("FDR analysis complete!\n")
