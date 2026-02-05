#!/usr/bin/env Rscript
################################################################################
# DEG_QQ.R
# Purpose: Generate QQ plots for DEG and DEG_SVs results
# Date: 02/05/2026
################################################################################

library(ggplot2)
library(dplyr)

# Configuration
deg_dir <- "DEG"
deg_svs_dir <- "DEG_SVs"
output_dir <- "QQ_plots"
dir.create(output_dir, showWarnings = FALSE)
today <- format(Sys.time(), "%m%d%Y")

################################################################################
# Part 1: DEG QQ plots (timepoints overlaid)
################################################################################

cat("Generating DEG QQ plots...\n")

contrasts <- c("BPA", "CRL", "MBP")
timepoints <- c(6, 24)

for (contrast in contrasts) {
  cat("Processing contrast:", contrast, "\n")
  
  # Collect data for both timepoints
  qq_data_list <- list()
  
  for (tp in timepoints) {
    # Find the file
    pattern <- paste0("_", contrast, "_T", tp, ".txt")
    file <- list.files(deg_dir, pattern = pattern, full.names = TRUE)
    
    if (length(file) == 0) {
      cat("Warning: No file found for", contrast, "T", tp, "\n")
      next
    }
    
    # Read results
    res <- read.table(file[1], header = TRUE, sep = "\t", row.names = 1)
    
    # Extract p-values and calculate observed/expected
    pvals <- res$P.Value
    pvals <- pvals[!is.na(pvals)]
    n <- length(pvals)
    
    observed <- -log10(sort(pvals))
    expected <- -log10(ppoints(n))
    
    qq_data_list[[paste0("T", tp)]] <- data.frame(
      expected = expected,
      observed = observed,
      timepoint = paste0("T", tp)
    )
  }
  
  # Combine data
  qq_data <- do.call(rbind, qq_data_list)
  
  # Create QQ plot
  p <- ggplot(qq_data, aes(x = expected, y = observed, color = timepoint)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
      title = paste("QQ Plot -", contrast),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_brewer(palette = "Set1")
  
  # Save plot
  filename <- file.path(output_dir, paste0(today, "_QQ_DEG_", contrast, ".png"))
  ggsave(filename, p, width = 8, height = 8, dpi = 300)
  cat("Saved:", filename, "\n")
}

################################################################################
# Part 2: DEG_SVs QQ plots (SV iterations overlaid)
################################################################################

cat("\nGenerating DEG_SVs QQ plots...\n")

max_svs <- 25

for (contrast in contrasts) {
  for (tp in timepoints) {
    cat("Processing contrast:", contrast, "T", tp, "\n")
    
    # Collect data for all SV iterations
    qq_data_list <- list()
    
    for (n_sv in 0:max_svs) {
      if (n_sv == 0) {
        sv_label <- "SV0"
      } else {
        sv_label <- paste0("SV1-", n_sv)
      }
      
      # Find the file
      pattern <- paste0("_", contrast, "_T", tp, "_", sv_label, ".txt")
      file <- list.files(deg_svs_dir, pattern = pattern, full.names = TRUE)
      
      if (length(file) == 0) {
        cat("Warning: No file found for", contrast, "T", tp, sv_label, "\n")
        next
      }
      
      # Read results
      res <- read.table(file[1], header = TRUE, sep = "\t", row.names = 1)
      
      # Extract p-values and calculate observed/expected
      pvals <- res$P.Value
      pvals <- pvals[!is.na(pvals)]
      n <- length(pvals)
      
      observed <- -log10(sort(pvals))
      expected <- -log10(ppoints(n))
      
      qq_data_list[[sv_label]] <- data.frame(
        expected = expected,
        observed = observed,
        sv_label = sv_label,
        n_sv = n_sv
      )
    }
    
    # Combine data
    qq_data <- do.call(rbind, qq_data_list)
    
    # Create QQ plot with gradient color
    p <- ggplot(qq_data, aes(x = expected, y = observed, color = n_sv, group = sv_label)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = paste("QQ Plot -", contrast, "T", tp),
        subtitle = "SV iterations (0-25)",
        x = "Expected -log10(p)",
        y = "Observed -log10(p)",
        color = "Number of SVs"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right"
      ) +
      scale_color_viridis_c(option = "turbo")
    
    # Save plot
    filename <- file.path(output_dir, paste0(today, "_QQ_DEG_SVs_", contrast, "_T", tp, ".png"))
    ggsave(filename, p, width = 10, height = 8, dpi = 300)
    cat("Saved:", filename, "\n")
  }
}

cat("\n=== Complete ===\n")
cat("Total plots created:", length(contrasts) + (length(contrasts) * length(timepoints)), "\n")
