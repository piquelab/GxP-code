#!/usr/bin/env Rscript
################################################################################
# PCA_vis.R
# Purpose: PCA visualization for filtered data
# Date: 02/05/2026
################################################################################

library(ggplot2)
library(gridExtra)
library(viridis)

### Configuration
input_file <- "GxP_Filtered_02052026.RData"
output_dir <- "PCA_analysis"
dir.create(output_dir, showWarnings = FALSE)
today <- format(Sys.time(), "%m%d%Y")

color_variables <- c("Timepoint", "Sequencer", "Batch_tx", "trimmed_dClean.dFastq", "sex", "Treatment")

### Load data
cat("Loading data...\n")
load(input_file)

### Combine both timepoints
cat("Combining timepoints...\n")

# Get voom data from both timepoints
voom_tp6 <- filtered_data$tp6$voom
voom_tp24 <- filtered_data$tp24$voom

# Find common genes
common_genes <- intersect(rownames(voom_tp6$E), rownames(voom_tp24$E))
cat("Common genes:", length(common_genes), "\n")

# Combine expression matrices
combined_expr <- cbind(
  voom_tp6$E[common_genes, ],
  voom_tp24$E[common_genes, ]
)

# Combine metadata
combined_metadata <- rbind(
  filtered_data$tp6$dge$samples,
  filtered_data$tp24$dge$samples
)

### Perform PCA
cat("Performing PCA...\n")
pca_result <- prcomp(t(combined_expr), center = TRUE, scale. = FALSE)

### Plotting function
plot_pca <- function(pca_result, metadata, color_var, pc_x, pc_y) {
  
  # Extract PC scores
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3],
    PC4 = pca_result$x[, 4],
    SampleID = rownames(pca_result$x)
  )
  
  # Merge with metadata
  pca_data <- merge(pca_data, metadata, by.x = "SampleID", by.y = "row.names")
  
  # Convert to factor if categorical
  if (color_var %in% c("Timepoint", "Sequencer", "Batch_tx", "sex", "Treatment")) {
    pca_data[[color_var]] <- as.factor(pca_data[[color_var]])
  }
  
  # Variance explained
  var_explained <- summary(pca_result)$importance[2, ]
  pc_x_var <- round(var_explained[pc_x] * 100, 1)
  pc_y_var <- round(var_explained[pc_y] * 100, 1)
  
  # Get PC columns
  pc_x_col <- paste0("PC", pc_x)
  pc_y_col <- paste0("PC", pc_y)
  
  # Create plot
  p <- ggplot(pca_data, aes(x = .data[[pc_x_col]], 
                            y = .data[[pc_y_col]], 
                            color = .data[[color_var]])) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = paste0("PC", pc_x, " vs PC", pc_y),
      subtitle = paste0("Both Timepoints | Colored by: ", color_var),
      x = paste0("PC", pc_x, " (", pc_x_var, "%)"),
      y = paste0("PC", pc_y, " (", pc_y_var, "%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  
  # Color scale based on variable type
  if (color_var == "trimmed_dClean.dFastq") {
    p <- p + scale_color_viridis_c()
  } else if (color_var == "Batch_tx") {
    # Use viridis discrete for many categories
    p <- p + scale_color_viridis_d(option = "turbo")
  } else {
    p <- p + scale_color_brewer(palette = "Set2")
  }
  
  return(p)
}

### Generate plots
cat("Generating plots...\n")

for (color_var in color_variables) {
  
  # Create PC1-2 plot
  p1 <- plot_pca(pca_result, combined_metadata, color_var, 1, 2)
  
  # Create PC3-4 plot
  p2 <- plot_pca(pca_result, combined_metadata, color_var, 3, 4)
  
  # Combine plots
  combined <- grid.arrange(p1, p2, nrow = 2)
  
  # Save combined plot
  filename <- file.path(
    output_dir,
    paste0(today, "_PCA_PC1-4_BothTimepoints_", color_var, ".png")
  )
  
  ggsave(filename, combined, width = 10, height = 14, dpi = 300)
}

### Variance explained plot (PC1-9)
cat("Generating variance explained plot...\n")

var_exp <- summary(pca_result)$importance[2, 1:9] * 100

var_data <- data.frame(
  PC = paste0("PC", 1:9),
  Variance = var_exp
)

p_var <- ggplot(var_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Variance Explained by Principal Components",
    subtitle = "Both Timepoints Combined",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_bw()

ggsave(file.path(output_dir, paste0(today, "_PCA_VarianceExplained.png")), 
       p_var, width = 10, height = 6, dpi = 300)

cat("\nDone! Plots saved to:", output_dir, "\n")
cat("Total plots:", length(color_variables) + 1, "\n")
