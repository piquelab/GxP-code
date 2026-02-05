#!/usr/bin/env Rscript
################################################################################
# Outlier_removal.R
# Purpose: Remove specified outlier samples from GxP RNA-seq dataset
# Date: 2025-01-09
################################################################################

# Input file
input_file <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/GxP_774Samples_20250803.RData"
load(input_file)

# Samples to remove - final
samples_to_remove <- c(
    "GxP12_02",
    "GxP12_03",
    "GxP12_05",
    "GxP12_06",
    "GxP37_02",
    "GxP37_05",
    "GxP11_08",
    "GxP11_11",
    "GxP22_07",
    "GxP22_09",
    "GxP22_10",
    "GxP22_12",
    "GxP24_07",
    "GxP24_10",
    "GxP21_20",
    "GxP21_23"
)

# Output files
today <- "01182026"
output_rdata <- paste0("GxP_SamplesRemoved_", today, ".RData")
output_removed_samples <- paste0("samples_removed_", today, ".csv")

# Check that required objects exist
required_objects <- c("c2", "m2", "m2.filt2", "m2.filt3")
setdiff(required_objects, ls()) # 0 = no missing objects

cat("Total samples in c2:", nrow(c2), "\n")
cat("Total genes in m2:", nrow(m2), "\n")

all(colnames(m2) == rownames(c2)) # Verify data alignment
samples_to_remove %in% rownames(c2) # Verify samples are in data

# Save removed sample info for posterity
removed_samples_data <- c2[samples_to_remove, ]
write.csv(removed_samples_data, 
          file = output_removed_samples,
          row.names = TRUE)
cat("Saved removed sample metadata to:", output_removed_samples, "\n")

# Remove from metadata
c2_original <- c2
c2 <- c2[!rownames(c2) %in% samples_to_remove, ]

# Remove from count matrices
m2_original <- m2
m2 <- m2[, !colnames(m2) %in% samples_to_remove]

m2.filt2_original <- m2.filt2
m2.filt2 <- m2.filt2[, !colnames(m2.filt2) %in% samples_to_remove]

m2.filt3_original <- m2.filt3
m2.filt3 <- m2.filt3[, !colnames(m2.filt3) %in% samples_to_remove]

cat("Samples remaining:", nrow(c2), "\n")
cat("Samples removed:", nrow(c2_original) - nrow(c2), "\n")

all(colnames(m2) == rownames(c2)) # Verify data alignment
all(colnames(m2.filt2) == rownames(c2)) # Verify data alignment
all(colnames(m2.filt3) == rownames(c2)) # Verify data alignment

save(c2, m2, m2.filt2, m2.filt3,
     file = output_rdata)

cat("\n=== SCRIPT COMPLETED SUCCESSFULLY ===\n")
cat("Files created:\n")
cat("  1.", output_rdata, "\n")
cat("  2.", output_removed_samples, "\n")
