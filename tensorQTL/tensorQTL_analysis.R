#!/usr/bin/env Rscript
# TensorQTL Results Analysis - FIXED VERSION

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(viridis)
library(corrplot)
library(pheatmap)  # Added for alternative heatmap

# =============================================================================
# CONFIGURATION
# =============================================================================
base_dir <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor"
input_dir <- file.path(base_dir, "tensorqtl_output_cis")
output_dir <- file.path(base_dir, "analysis_results")
figures_dir <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Define treatment conditions
conditions <- c("BPA_100nM_24", "BPA_100nM_6", "EtOH_24", "EtOH_6",
                "H2O_24", "H2O_6", "MBP_500nM_24", "MBP_500nM_6")

cat("Starting tensorQTL analysis for", length(conditions), "conditions\n")

# =============================================================================
# 1. LOAD AND PROCESS DATA
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
    qtl_data[[condition]] <- qtl_results

    # Calculate summary statistics
    n_tests <- nrow(qtl_results)
    n_genes <- length(unique(qtl_results$phenotype_id))
    n_variants <- length(unique(qtl_results$variant_id))

    # Add FDR correction
    qtl_results$pval_fdr <- p.adjust(qtl_results$pval_nominal, method = "fdr")

    # Count significant eQTLs using standard R filtering (not data.table syntax)
    sig_threshold <- 0.05 / n_tests
    sig_mask_bonf <- qtl_results$pval_nominal < sig_threshold & !is.na(qtl_results$pval_nominal)
    sig_mask_fdr <- qtl_results$pval_fdr < 0.05 & !is.na(qtl_results$pval_fdr)
    sig_mask_001 <- qtl_results$pval_nominal < 0.001 & !is.na(qtl_results$pval_nominal)
    sig_mask_0001 <- qtl_results$pval_nominal < 0.0001 & !is.na(qtl_results$pval_nominal)

    n_sig_bonferroni <- sum(sig_mask_bonf)
    n_sig_fdr <- sum(sig_mask_fdr)
    n_sig_pval001 <- sum(sig_mask_001)
    n_sig_pval0001 <- sum(sig_mask_0001)

    # Count eGenes
    egenes_bonferroni <- length(unique(qtl_results$phenotype_id[sig_mask_bonf]))
    egenes_fdr <- length(unique(qtl_results$phenotype_id[sig_mask_fdr]))
    egenes_001 <- length(unique(qtl_results$phenotype_id[sig_mask_001]))
    egenes_0001 <- length(unique(qtl_results$phenotype_id[sig_mask_0001]))

    summary_row <- data.frame(
      condition = condition,
      n_tests = n_tests,
      n_genes = n_genes,
      n_variants = n_variants,
      egenes_bonferroni = egenes_bonferroni,
      egenes_fdr = egenes_fdr,
      egenes_p001 = egenes_001,
      egenes_p0001 = egenes_0001,
      eqtls_bonferroni = n_sig_bonferroni,
      eqtls_fdr = n_sig_fdr,
      eqtls_p001 = n_sig_pval001,
      eqtls_p0001 = n_sig_pval0001,
      bonferroni_threshold = sig_threshold
    )

    qtl_summary <- rbind(qtl_summary, summary_row)
  } else {
    cat("WARNING: File not found:", qtl_file, "\n")
  }
}

cat("Data loading complete. Total tests:", sum(qtl_summary$n_tests), "\n")

# =============================================================================
# 2. SUMMARY STATISTICS TABLE  
# =============================================================================
cat("Generating summary statistics...\n")

# Add treatment and timepoint columns
qtl_summary$treatment <- sapply(strsplit(qtl_summary$condition, "_"), function(x) {
  if(length(x) >= 1) x[1] else "unknown"
})

qtl_summary$timepoint <- sapply(strsplit(qtl_summary$condition, "_"), function(x) {
  if(length(x) >= 2) {
    paste(x[2:length(x)], collapse = "_")
  } else {
    "unknown"
  }
})

# Write summary table
fwrite(qtl_summary, file.path(output_dir, "eQTL_summary_by_condition.txt"),
       sep = '\t', quote = FALSE, row.names = FALSE)

# =============================================================================
# 3. VISUALIZATION: SUMMARY BARPLOTS  
# =============================================================================
cat("Creating summary visualizations...\n")

# eGenes by condition (FDR corrected)
p1 <- ggplot(qtl_summary, aes(x = condition, y = egenes_fdr)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = egenes_fdr), vjust = -0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of eGenes by Condition (FDR < 0.05)",
       x = "Condition", y = "Number of eGenes") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figures_dir, "eGenes_by_condition_FDR.png"), p1,
       width = 10, height = 6, dpi = 300)

# eGenes by condition (nominal p < 0.001 for comparison)
p1_nominal <- ggplot(qtl_summary, aes(x = condition, y = egenes_p001)) +
  geom_bar(stat = "identity", fill = "lightblue", alpha = 0.7) +
  geom_text(aes(label = egenes_p001), vjust = -0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of eGenes by Condition (nominal p < 0.001)",
       x = "Condition", y = "Number of eGenes") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figures_dir, "eGenes_by_condition_nominal.png"), p1_nominal,
       width = 10, height = 6, dpi = 300)

# eQTLs by condition (FDR corrected)
p2 <- ggplot(qtl_summary, aes(x = condition, y = eqtls_fdr)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  geom_text(aes(label = eqtls_fdr), vjust = -0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Significant eQTLs by Condition (FDR < 0.05)",
       x = "Condition", y = "Number of eQTLs") +   
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figures_dir, "eQTLs_by_condition_FDR.png"), p2,
       width = 10, height = 6, dpi = 300)

# eQTLs by condition (nominal p < 0.001 for comparison)
p2_nominal <- ggplot(qtl_summary, aes(x = condition, y = eqtls_p001)) +
  geom_bar(stat = "identity", fill = "lightgreen", alpha = 0.7) +
  geom_text(aes(label = eqtls_p001), vjust = -0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Significant eQTLs by Condition (nominal p < 0.001)",
       x = "Condition", y = "Number of eQTLs") +   
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figures_dir, "eQTLs_by_condition_nominal.png"), p2_nominal,
       width = 10, height = 6, dpi = 300)

# =============================================================================
# 4. QQ PLOTS FOR EACH CONDITION
# =============================================================================
cat("Creating QQ plots...\n")

create_qq_plot <- function(condition_data, condition_name, output_path) {
  pvals <- condition_data$pval_nominal[!is.na(condition_data$pval_nominal)]
  pvals <- pvals[pvals > 0]
  pvals <- sort(pvals)

  n <- length(pvals)
  expected <- (1:n) / n
  pvals[pvals < 1e-20] <- 1e-20

  qq_data <- data.frame(
    expected = -log10(expected),
    observed = -log10(pvals)
  )

  ci <- 0.95
  qq_data$clower <- -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n - (1:n) + 1))
  qq_data$cupper <- -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n - (1:n) + 1))

  p <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2, fill = "gray") +
    geom_point(alpha = 0.6, size = 0.8) +  
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_classic() +
    labs(title = paste("QQ Plot -", condition_name),
         x = expression(Expected -log[10](p)),
         y = expression(Observed -log[10](p))) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))

  ggsave(output_path, p, width = 8, height = 8, dpi = 300) 
  return(p)
}

for (condition in conditions) {
  if (condition %in% names(qtl_data)) {  
    output_file <- file.path(figures_dir, paste0("QQ_plot_", condition, ".png"))
    create_qq_plot(qtl_data[[condition]], condition, output_file)
  }
}

# =============================================================================
# 5. TREATMENT COMPARISON ANALYSIS (FIXED WITH STANDARD R SYNTAX)
# =============================================================================
cat("Performing treatment comparisons...\n")

# Extract significant eQTLs using FDR correction (primary) and nominal p-value (secondary)
sig_eqtls_fdr <- list()
sig_eqtls_nominal <- list()

for (condition in conditions) {
  if (condition %in% names(qtl_data)) {
    cat("Processing overlaps for:", condition, "\n")

    # Use standard R filtering instead of data.table syntax
    current_data <- qtl_data[[condition]]
    
    # FDR-corrected significant eQTLs (primary analysis)
    sig_mask_fdr <- current_data$pval_fdr < 0.05 & !is.na(current_data$pval_fdr)
    sig_data_fdr <- current_data[sig_mask_fdr, ]
    
    # Nominal p-value significant eQTLs (for comparison)
    sig_mask_nominal <- current_data$pval_nominal < 0.001 & !is.na(current_data$pval_nominal)
    sig_data_nominal <- current_data[sig_mask_nominal, ]

    if (nrow(sig_data_fdr) > 0) {
      sig_eqtls_fdr[[condition]] <- unique(sig_data_fdr$phenotype_id)
    } else {
      sig_eqtls_fdr[[condition]] <- character(0)
    }
    
    if (nrow(sig_data_nominal) > 0) {
      sig_eqtls_nominal[[condition]] <- unique(sig_data_nominal$phenotype_id)
    } else {
      sig_eqtls_nominal[[condition]] <- character(0)
    }
    
    cat("  Found", length(sig_eqtls_fdr[[condition]]), "eGenes (FDR < 0.05)\n")
    cat("  Found", length(sig_eqtls_nominal[[condition]]), "eGenes (nominal p < 0.001)\n")
  }
}

# Create overlap matrices for both FDR and nominal thresholds
valid_conditions <- names(sig_eqtls_fdr)

# FDR-based overlap matrix (primary analysis)
overlap_matrix_fdr <- matrix(0, nrow = length(valid_conditions), ncol = length(valid_conditions))
rownames(overlap_matrix_fdr) <- valid_conditions
colnames(overlap_matrix_fdr) <- valid_conditions

# Nominal-based overlap matrix (for comparison)
overlap_matrix_nominal <- matrix(0, nrow = length(valid_conditions), ncol = length(valid_conditions))
rownames(overlap_matrix_nominal) <- valid_conditions
colnames(overlap_matrix_nominal) <- valid_conditions

cat("Building overlap matrices...\n")

# Build FDR-based overlap matrix (primary)
for (i in 1:length(valid_conditions)) {
  for (j in 1:length(valid_conditions)) {
    cond_i <- valid_conditions[i]
    cond_j <- valid_conditions[j]
    overlap_matrix_fdr[i, j] <- length(intersect(sig_eqtls_fdr[[cond_i]], sig_eqtls_fdr[[cond_j]]))
  }
}

# Build nominal-based overlap matrix (for comparison)
for (i in 1:length(valid_conditions)) {
  for (j in 1:length(valid_conditions)) {
    cond_i <- valid_conditions[i]
    cond_j <- valid_conditions[j]
    overlap_matrix_nominal[i, j] <- length(intersect(sig_eqtls_nominal[[cond_i]], sig_eqtls_nominal[[cond_j]]))
  }
}

# Save both overlap matrices
write.table(overlap_matrix_fdr, file.path(output_dir, "eGene_overlap_matrix_FDR.txt"),
            sep = '\t', quote = FALSE)
write.table(overlap_matrix_nominal, file.path(output_dir, "eGene_overlap_matrix_nominal.txt"),
            sep = '\t', quote = FALSE)

# FIXED: Convert overlap matrices to correlation/similarity matrices for corrplot
if (length(valid_conditions) > 1) {
  
  # === FDR-BASED ANALYSIS (PRIMARY) ===
  # Check if we have any FDR-significant results
  total_fdr_egenes <- sum(sapply(sig_eqtls_fdr, length))
  
  if (total_fdr_egenes > 0) {
    cat("Creating FDR-based heatmaps...\n")
    
    # Method 1: Convert to Jaccard similarity (values between 0 and 1)
    jaccard_matrix_fdr <- matrix(0, nrow = nrow(overlap_matrix_fdr), ncol = ncol(overlap_matrix_fdr))
    rownames(jaccard_matrix_fdr) <- rownames(overlap_matrix_fdr)
    colnames(jaccard_matrix_fdr) <- colnames(overlap_matrix_fdr)
    
    for (i in 1:nrow(overlap_matrix_fdr)) {
      for (j in 1:ncol(overlap_matrix_fdr)) {
        cond_i <- valid_conditions[i]
        cond_j <- valid_conditions[j]
        intersection <- overlap_matrix_fdr[i, j]
        union_size <- length(sig_eqtls_fdr[[cond_i]]) + length(sig_eqtls_fdr[[cond_j]]) - intersection
        jaccard_matrix_fdr[i, j] <- if (union_size > 0) intersection / union_size else 0
      }
    }
    
    # Create corrplot with FDR-based Jaccard similarity
    png(file.path(figures_dir, "eGene_jaccard_similarity_heatmap_FDR.png"),
        width = 12, height = 10, units = "in", res = 300)
    corrplot(jaccard_matrix_fdr, method = "color", type = "full",
             order = "hclust", tl.cex = 0.8, tl.col = "black",
             addCoef.col = "black", number.cex = 0.7,
             title = "eGene Jaccard Similarity Between Conditions (FDR < 0.05)", 
             mar = c(0,0,2,0), col = viridis(100))
    dev.off()
    
    # Create FDR-based raw count heatmap
    png(file.path(figures_dir, "eGene_overlap_counts_heatmap_FDR.png"),
        width = 12, height = 10, units = "in", res = 300)
    pheatmap(overlap_matrix_fdr, 
             display_numbers = TRUE, 
             number_color = "white",
             color = viridis(100),
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             main = "eGene Overlap Counts Between Conditions (FDR < 0.05)")
    dev.off()
  } else {
    cat("WARNING: No FDR-significant eGenes found across all conditions.\n")
    cat("Skipping FDR-based heatmaps. Consider using a less stringent FDR threshold.\n")
    
    # Create a note file explaining the issue
    warning_msg <- paste0(
      "FDR Analysis Warning\n",
      "===================\n\n",
      "No eGenes were significant at FDR < 0.05 across all conditions.\n",
      "This suggests that the FDR correction is very stringent for your data.\n\n",
      "Recommendations:\n",
      "1. Use nominal p-value results (p < 0.001) for exploratory analysis\n",
      "2. Consider a less stringent FDR threshold (e.g., FDR < 0.10 or 0.20)\n",
      "3. Review your sample size and effect sizes\n",
      "4. Consider using tensorQTL's own q-values if available in your output files\n\n",
      "The analysis will continue with nominal p-value thresholds.\n"
    )
    writeLines(warning_msg, file.path(output_dir, "FDR_WARNING.txt"))
  }
  
  # === NOMINAL-BASED ANALYSIS (FOR COMPARISON) ===
  # Convert nominal overlap matrix to Jaccard similarity
  jaccard_matrix_nominal <- matrix(0, nrow = nrow(overlap_matrix_nominal), ncol = ncol(overlap_matrix_nominal))
  rownames(jaccard_matrix_nominal) <- rownames(overlap_matrix_nominal)
  colnames(jaccard_matrix_nominal) <- colnames(overlap_matrix_nominal)
  
  for (i in 1:nrow(overlap_matrix_nominal)) {
    for (j in 1:ncol(overlap_matrix_nominal)) {
      cond_i <- valid_conditions[i]
      cond_j <- valid_conditions[j]
      intersection <- overlap_matrix_nominal[i, j]
      union_size <- length(sig_eqtls_nominal[[cond_i]]) + length(sig_eqtls_nominal[[cond_j]]) - intersection
      jaccard_matrix_nominal[i, j] <- if (union_size > 0) intersection / union_size else 0
    }
  }
  
  # Create corrplot with nominal-based Jaccard similarity
  png(file.path(figures_dir, "eGene_jaccard_similarity_heatmap_nominal.png"),
      width = 12, height = 10, units = "in", res = 300)
  corrplot(jaccard_matrix_nominal, method = "color", type = "full",
           order = "hclust", tl.cex = 0.8, tl.col = "black",
           addCoef.col = "black", number.cex = 0.7,
           title = "eGene Jaccard Similarity Between Conditions (nominal p < 0.001)", 
           mar = c(0,0,2,0), col = viridis(100))
  dev.off()
  
  # Create nominal-based raw count heatmap
  png(file.path(figures_dir, "eGene_overlap_counts_heatmap_nominal.png"),
      width = 12, height = 10, units = "in", res = 300)
  pheatmap(overlap_matrix_nominal, 
           display_numbers = TRUE, 
           number_color = "white",
           color = viridis(100),
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           main = "eGene Overlap Counts Between Conditions (nominal p < 0.001)")
  dev.off()
  
  cat("Created both FDR-corrected and nominal threshold heatmaps\n")
} else {
  cat("Only one condition found, skipping overlap analysis\n")
}

# =============================================================================
# 6. SAVE SIGNIFICANT RESULTS (FIXED)
# =============================================================================
cat("Saving significant results...\n")

for (condition in conditions) {
  if (condition %in% names(qtl_data)) {
    cat("Saving results for:", condition, "\n")

    # Use standard R filtering 
    current_data <- qtl_data[[condition]]
    
    # Save FDR-corrected significant results (primary analysis)
    sig_mask_fdr <- current_data$pval_fdr < 0.05 & !is.na(current_data$pval_fdr)
    sig_results_fdr <- current_data[sig_mask_fdr, ]

    if (nrow(sig_results_fdr) > 0) {
      # Sort by FDR-corrected p-value
      sig_results_fdr <- sig_results_fdr[order(sig_results_fdr$pval_fdr), ]

      # Save full FDR results
      output_file_fdr <- file.path(output_dir, paste0(condition, "_significant_eQTLs_FDR005.txt"))
      fwrite(sig_results_fdr, output_file_fdr, sep = '\t', quote = FALSE, row.names = FALSE)

      # Save FDR gene-variant pairs using standard R syntax
      pairs_fdr <- sig_results_fdr[, c("phenotype_id", "variant_id")]
      pairs_file_fdr <- file.path(output_dir, paste0(condition, "_significant_pairs_FDR005.txt"))
      fwrite(pairs_fdr, pairs_file_fdr, sep = '\t', quote = FALSE, row.names = FALSE)

      cat("  Saved", nrow(sig_results_fdr), "significant eQTLs (FDR < 0.05)\n")
    } else {
      cat("  No significant FDR results for", condition, "\n")
    }
    
    # Also save nominal results for comparison
    sig_mask_nominal <- current_data$pval_nominal < 0.001 & !is.na(current_data$pval_nominal)
    sig_results_nominal <- current_data[sig_mask_nominal, ]

    if (nrow(sig_results_nominal) > 0) {
      # Sort by nominal p-value
      sig_results_nominal <- sig_results_nominal[order(sig_results_nominal$pval_nominal), ]

      # Save full nominal results
      output_file_nominal <- file.path(output_dir, paste0(condition, "_significant_eQTLs_nominal_p001.txt"))
      fwrite(sig_results_nominal, output_file_nominal, sep = '\t', quote = FALSE, row.names = FALSE)

      # Save nominal gene-variant pairs using standard R syntax  
      pairs_nominal <- sig_results_nominal[, c("phenotype_id", "variant_id")]
      pairs_file_nominal <- file.path(output_dir, paste0(condition, "_significant_pairs_nominal_p001.txt"))
      fwrite(pairs_nominal, pairs_file_nominal, sep = '\t', quote = FALSE, row.names = FALSE)

      cat("  Saved", nrow(sig_results_nominal), "significant eQTLs (nominal p < 0.001)\n")
    } else {
      cat("  No significant nominal results for", condition, "\n")
    }
  }
}

# =============================================================================
# 7. FINAL SUMMARY REPORT
# =============================================================================
cat("Generating final summary report...\n")

total_egenes_fdr <- sum(qtl_summary$egenes_fdr)
total_egenes_nominal <- sum(qtl_summary$egenes_p001)
total_eqtls_fdr <- sum(qtl_summary$eqtls_fdr)
total_eqtls_nominal <- sum(qtl_summary$eqtls_p001)
total_tests <- sum(qtl_summary$n_tests)

summary_report <- paste0(
  "TensorQTL Analysis Summary\n",
  paste(rep("=", 50), collapse = ""), "\n",
  "Total conditions analyzed: ", length(conditions), "\n",
  "Total statistical tests: ", format(total_tests, big.mark = ","), "\n\n",
  "=== FDR-CORRECTED RESULTS (PRIMARY) ===\n",
  "Total eGenes (FDR < 0.05): ", total_egenes_fdr, "\n",
  "Total significant eQTLs (FDR < 0.05): ", format(total_eqtls_fdr, big.mark = ","), "\n\n",
  "=== NOMINAL RESULTS (FOR COMPARISON) ===\n",
  "Total eGenes (nominal p < 0.001): ", total_egenes_nominal, "\n",
  "Total significant eQTLs (nominal p < 0.001): ", format(total_eqtls_nominal, big.mark = ","), "\n\n",
  "Results saved to: ", output_dir, "\n",
  "Figures saved to: ", figures_dir, "\n\n",
  "Files generated:\n",
  "- eQTL_summary_by_condition.txt: Summary statistics table\n",
  "- *_significant_eQTLs_FDR005.txt: FDR-corrected significant eQTL results per condition\n",
  "- *_significant_pairs_FDR005.txt: FDR-corrected gene-variant pairs per condition\n",
  "- *_significant_eQTLs_nominal_p001.txt: Nominal significant eQTL results per condition\n",
  "- *_significant_pairs_nominal_p001.txt: Nominal gene-variant pairs per condition\n",
  "- eGene_overlap_matrix_FDR.txt: FDR-corrected overlap counts between conditions\n",
  "- eGene_overlap_matrix_nominal.txt: Nominal overlap counts between conditions\n",
  "- figures/: QQ plots and summary visualizations (both FDR and nominal)\n"
)

writeLines(summary_report, file.path(output_dir, "ANALYSIS_SUMMARY.txt"))
cat(summary_report)

cat("Analysis complete!\n")
