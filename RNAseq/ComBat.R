#!/usr/bin/env Rscript
################################################################################
# ComBat.R
# Purpose: Re-run ComBat_seq batch correction on multiple batch variables after sample removal. Removes old ComBat data and saves clean file.
# Date: 2025-01-18
################################################################################

library(sva)

### Input file (output from outlier_removal.R)
input_file <- "GxP_SamplesRemoved_01182026.RData"
load(input_file)

### Output file
today <- "01182026"
output_file <- paste0("GxP_ComBatCorrected_", today, ".RData")

### Batch variables to correct for
batch_variables <- c("Sequencer", "Batch_tx")

### Corresponding output names for corrected matrices
output_names <- c("m2.seq", "m2.batch")

cat("\nData dimensions:\n")
cat("  Samples (c2):", nrow(c2), "\n")
cat("  Genes (m2):", nrow(m2), "\n")
cat("  Samples in counts (m2):", ncol(m2), "\n")

### Remove old ComBat-corrected matrices if they exist
rm(list = c("m2.filt2", "m2.filt3"), 
   envir = .GlobalEnv)

### Check batch distributions
cat("Batch Variable Distributions \n")
for (batch_var in batch_variables) {
  cat("\n", batch_var, ":\n", sep = "")
  print(table(c2[[batch_var]]))
  
  cat("  vs Treatment:\n")
  print(table(c2[[batch_var]], c2$Treatment))
}

############## COMBAT ####################

### Convert to matrix once
m2_matrix <- as.matrix(m2)

### Storage for results
combat_results <- list()

### Run ComBat for each batch variable
for (i in seq_along(batch_variables)) {
  batch_var <- batch_variables[i]
  output_name <- output_names[i]
  
  cat("\nRunning ComBat_seq for:", batch_var, "→", output_name, "\n")
  
  tryCatch({
    corrected_matrix <- ComBat_seq(
      counts = m2_matrix,
      batch = c2[[batch_var]],
      group = NULL,
      full_mod = TRUE
    )
    assign(output_name, corrected_matrix, envir = .GlobalEnv)
    cat("✓ Completed\n")
  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n")
  })
}


### Verification
for (name in output_names) {
  if (exists(name)) {
    mat <- get(name)
    aligned <- all(colnames(mat) == rownames(c2))
    corr <- cor(m2[,1], mat[,1])
    cat(sprintf("  %-15s %d×%d  aligned:%s  cor=%.3f\n", 
                name, nrow(mat), ncol(mat), 
                ifelse(aligned, "✓", "✗"), corr))
  }
}

### Save results
objects_to_save <- c("c2", "m2", output_names[sapply(output_names, exists)])
save(list = objects_to_save, file = output_file)
cat("✓ Saved to:", output_file, "\n")
