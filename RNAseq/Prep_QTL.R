#!/usr/bin/env Rscript
################################################################################
# Prep_QTL.R
# Purpose: Prepare expression data for QTL mapping with SVA
# Date: 02/05/2026
################################################################################

library(edgeR)
library(limma)
library(preprocessCore)
library(annotables)
library(sva)
library(dplyr)

# Configuration
BASE <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026"
INPUT_RDATA <- file.path(BASE, "GxP_Filtered_02052026.RData")
VCF_FILE <- file.path(BASE, "GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz")
OUTPUT_DIR <- BASE
SV_DIR <- file.path(BASE, "Norm_SVs")
dir.create(SV_DIR, showWarnings = FALSE)

TIMEPOINTS <- c(6, 24)
PAIRS <- list(
  list(treatments = c("BPA_100nM", "H2O"), name_suffix = c("BPA", "H2O"), pair_name = "BPA_H2O"),
  list(treatments = c("MBP_500nM", "EtOH"), name_suffix = c("MBP", "EtOH"), pair_name = "MBP_EtOH")
)

# Load data
cat("Loading data...\n")
load(INPUT_RDATA)

# Load annotation
anno <- grch38 %>% 
  dplyr::select(ensgene, chr, start, end, biotype) %>% 
  dplyr::filter(biotype == "protein_coding", chr %in% c(1:22)) %>% 
  unique()
anno <- anno[!duplicated(anno$ensgene), ]

# Get genotyped samples
cat("Reading VCF samples...\n")
vcf_header <- system(paste("zcat", VCF_FILE, "| grep -m1 '^#CHROM'"), intern = TRUE)
vcf_samples <- unlist(strsplit(vcf_header, "\t"))[-c(1:9)]

# Process each timepoint
for (tp in TIMEPOINTS) {
  cat("\n=== Timepoint", tp, "===\n")
  
  # Get data for this timepoint
  tp_data <- filtered_data[[paste0("tp", tp)]]
  dge <- tp_data$dge
  samples <- dge$samples
  
  # Voom normalize
  cat("Voom normalizing...\n")
  design_basic <- model.matrix(~ Treatment + dbGaP_ID + trimmed_dClean.dFastq, data = samples)
  voom_data <- voom(dge, design = design_basic, plot = FALSE)
  
  # Regress out trimmed_dClean.dFastq
  cat("Regressing out trimmed_dClean.dFastq...\n")
  design_regress <- model.matrix(~ trimmed_dClean.dFastq, data = samples)
  fit_regress <- lmFit(voom_data, design_regress)
  residuals <- residuals(fit_regress, voom_data)
  
  ############################################################################
  # Process each treatment pair
  ############################################################################
  
  for (pair in PAIRS) {
    cat("\n--- Treatment pair:", pair$pair_name, "---\n")
    
    # Get samples for this pair
    pair_samples <- samples[samples$Treatment %in% pair$treatments, ]
    pair_expr <- residuals[, rownames(pair_samples)]
    
    # Quantile normalize
    cat("Quantile normalizing...\n")
    qnorm_expr <- normalize.quantiles(pair_expr)
    rownames(qnorm_expr) <- rownames(pair_expr)
    colnames(qnorm_expr) <- colnames(pair_expr)
    
    # Run SVA
    cat("Running SVA...\n")
    pair_samples$Treatment <- factor(pair_samples$Treatment, levels = pair$treatments)
    pair_samples$dbGaP_ID <- make.names(pair_samples$dbGaP_ID)
    
    # Full and null models
    mod <- model.matrix(~ Treatment + dbGaP_ID, data = pair_samples)
    mod0 <- model.matrix(~ dbGaP_ID, data = pair_samples)
    
    # Estimate number of SVs
    n.sv <- num.sv(qnorm_expr, mod, method = "be")
    cat("Estimated SVs:", n.sv, "\n")
    
    # Run SVA
    svobj <- sva(qnorm_expr, mod, mod0, n.sv = n.sv)
    cat("Identified", svobj$n.sv, "surrogate variables\n")
    
    # SV matrix (with RNA sample IDs)
    sv_matrix <- svobj$sv
    rownames(sv_matrix) <- colnames(qnorm_expr)
    colnames(sv_matrix) <- paste0("SV", 1:svobj$n.sv)
    
    ##########################################################################
    # Split SVs by condition and save with dbGaP_IDs
    ##########################################################################
    
    for (i in seq_along(pair$treatments)) {
      treatment <- pair$treatments[i]
      condition_name <- pair$name_suffix[i]
      
      cat("\nSaving SVs for condition:", condition_name, "\n")
      
      # Get samples for this condition
      cond_samples <- pair_samples[pair_samples$Treatment == treatment, ]
      
      # Map to dbGaP IDs
      sample_mapping <- data.frame(
        rna_id = rownames(cond_samples),
        dbgap_id = cond_samples$dbGaP_ID
      )
      
      # Convert periods to hyphens to match VCF format
      sample_mapping$dbgap_id_vcf <- gsub("\\.", "-", sample_mapping$dbgap_id)
      
      # Filter to genotyped samples
      genotyped <- sample_mapping$dbgap_id_vcf %in% vcf_samples
      sample_mapping_filt <- sample_mapping[genotyped, ]
      
      if (nrow(sample_mapping_filt) == 0) {
        warning("No genotyped samples for ", condition_name, "_T", tp)
        next
      }
      
      cat("Genotyped samples:", nrow(sample_mapping_filt), "/", nrow(sample_mapping), "\n")
      
      # Subset SV matrix to this condition's samples
      sv_cond <- sv_matrix[sample_mapping_filt$rna_id, , drop = FALSE]
      
      # Replace rownames with dbGaP IDs (VCF format)
      rownames(sv_cond) <- sample_mapping_filt$dbgap_id_vcf
      
      # Transpose for tensorQTL format (SVs as rows, samples as columns)
      sv_cond_t <- t(sv_cond)
      
      # Add id column
      sv_out <- data.frame(id = rownames(sv_cond_t), sv_cond_t, check.names = FALSE)
      
      # Save condition-specific SV file
      sv_file <- file.path(SV_DIR, paste0("SV_", condition_name, "_T", tp, ".txt"))
      write.table(sv_out, sv_file, sep = "\t", quote = FALSE, row.names = FALSE)
      cat("Saved:", sv_file, "\n")
      cat("Dimensions:", nrow(sv_out), "SVs x", ncol(sv_out) - 1, "samples\n")
    }
    
    ##########################################################################
    # Split into conditions and write BED files
    ##########################################################################
    
    for (i in seq_along(pair$treatments)) {
      treatment <- pair$treatments[i]
      condition_name <- pair$name_suffix[i]
      
      cat("\nProcessing BED for condition:", condition_name, "\n")
      
      # Get samples for this condition
      cond_samples <- pair_samples[pair_samples$Treatment == treatment, ]
      cond_expr <- qnorm_expr[, rownames(cond_samples)]
      
      # Map to dbGaP IDs and filter to genotyped
      sample_mapping <- data.frame(
        rna_id = rownames(cond_samples),
        dbgap_id = cond_samples$dbGaP_ID
      )
      
      # Convert periods to hyphens to match VCF format
      sample_mapping$dbgap_id_vcf <- gsub("\\.", "-", sample_mapping$dbgap_id)
      
      genotyped <- sample_mapping$dbgap_id_vcf %in% vcf_samples
      sample_mapping_filt <- sample_mapping[genotyped, ]
      
      if (nrow(sample_mapping_filt) == 0) {
        warning("No genotyped samples for ", condition_name, "_T", tp)
        next
      }
      
      # Filter expression and rename to dbGaP IDs (VCF format)
      cond_expr_filt <- cond_expr[, sample_mapping_filt$rna_id, drop = FALSE]
      colnames(cond_expr_filt) <- sample_mapping_filt$dbgap_id_vcf
      
      # Add gene coordinates
      gene_data <- anno[match(rownames(cond_expr_filt), anno$ensgene), ]
      
      # Create BED: #chr start end gene_id sample1 sample2 ...
      bed_data <- data.frame(
        chr = gene_data$chr,
        start = gene_data$start - 1,  # 0-based
        end = gene_data$start,
        gene_id = rownames(cond_expr_filt),
        cond_expr_filt,
        check.names = FALSE
      )
      
      # Sort by position
      bed_data <- bed_data[order(as.numeric(bed_data$chr), bed_data$start), ]
      
      # Write BED
      output_file <- file.path(OUTPUT_DIR, paste0("GxP-eQTL_", condition_name, "_T", tp, "_qnorm.bed.gz"))
      gz_con <- gzfile(output_file, "w")
      writeLines(paste0("#", paste(colnames(bed_data), collapse = "\t")), gz_con)
      write.table(bed_data, gz_con, sep = "\t", quote = FALSE, 
                  row.names = FALSE, col.names = FALSE)
      close(gz_con)
      
      cat("Saved BED:", output_file, "\n")
      cat("Genes:", nrow(bed_data), "Samples:", ncol(cond_expr_filt), "\n")
    }
  }
}

cat("\n=== Complete ===\n")
output_files <- list.files(OUTPUT_DIR, pattern = "^GxP-eQTL_.*_qnorm\\.bed\\.gz$")
cat("Generated", length(output_files), "BED files\n")

sv_files <- list.files(SV_DIR, pattern = "^SV_.*\\.txt$")
cat("Generated", length(sv_files), "condition-specific SV files\n")
