#!/usr/bin/env Rscript
################################################################################
# Filter_expression.R
# Purpose: Separate by time point, filter for expressed genes, TMM normalize, filter for autosomal coding genes, voom normalize, regress out covariate, split by pair, quantile normalize, inverse normal transform, calculate SVs.
################################################################################

library(edgeR)
library(preprocessCore)
library(annotables)
library(dplyr)
library(limma)

# Configure inputs
input   <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/starter_files/GxP_SamplesRemoved_01182026.RData"
output  <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/GxP_Filtered.RData"
VCF_FILE <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/starter_files/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz"
OUTPUT_DIR <- "/rs/rs_grp_gxp/RNAseq_analysis/phthalates/starter_files"
PC_DIR     <- file.path(OUTPUT_DIR, "PCs")
dir.create(PC_DIR, showWarnings = FALSE)

TIMEPOINTS  <- c(6, 24)
TREATMENTS  <- c("MBP_500nM", "EtOH")
# Short labels for file naming: index matches TREATMENTS
TX_LABELS   <- c("MBP", "EtOH")

# Load data and annotations
load(input)

# Read VCF sample IDs
cat("Reading VCF samples...\n")
vcf_header  <- system(paste("zcat", VCF_FILE, "| grep -m1 '^#CHROM'"), intern = TRUE)
vcf_samples <- unlist(strsplit(vcf_header, "\t"))[-c(1:9)]

# Prime storage for results
filtered_data <- list()

  # Filter and clean gene IDs
  colData <- c2[c2$Treatment %in% TREATMENTS, ]
  counts  <- m2[, rownames(colData)]
  rownames(counts) <- sub("\\..*", "", rownames(counts))

  # Build DGE, filter, normalize
  dge <- DGEList(counts = counts, samples = colData)
  
  # CPM-based filtering (GTEx-style: CPM >= 0.1 in >= 20% of samples)
  cpm_mat <- cpm(dge)
  keep <- rowSums(cpm_mat >= 0.1) >= (0.2 * ncol(dge))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  
  # TMM normalization
  dge  <- calcNormFactors(dge, method = "TMM")
  
  # Filter for protein-coding autosomal genes
  anno <- grch38 %>%
    dplyr::select(ensgene, chr, start, end, strand, biotype) %>%
    dplyr::filter(biotype == "protein_coding", chr %in% c(1:22)) %>%
    unique()
  anno <- anno[!duplicated(anno$ensgene), ]
  
  dge  <- dge[rownames(dge) %in% anno$ensgene, ]
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Set treatment levels
  dge$samples$Batch_tx   <- factor(dge$samples$Batch_tx)
  dge$samples$Sequencer  <- factor(dge$samples$Sequencer)
  dge$samples$Treatment <- factor(dge$samples$Treatment, levels = TREATMENTS)
  dge$samples$dbGaP_ID   <- make.names(dge$samples$dbGaP_ID)

  # SAVE FILTERED DGE FOR DEG ANALYSIS - voom run per timepoint in Limma.R
  filtered_data <- list(dge = dge)
  save(filtered_data, anno, c2, file = output)
  cat("Saved RData to:", output, "\n")

  # Regress out sequencing depth covariate
  cat("Regressing out trimmed_dClean.dFastq...\n")
  logcpm      <- cpm(dge, log = TRUE, prior.count = 0.5)
  fit_regress <- lmFit(logcpm, model.matrix(~ trimmed_dClean.dFastq, data = dge$samples))
  residuals   <- residuals(fit_regress, logcpm)

 # Per-condition: subset by treatment x timepoint -> INT -> SVA -> save
for (tp in TIMEPOINTS) {
  for (tx_idx in seq_along(TREATMENTS)) {
 
    treatment      <- TREATMENTS[tx_idx]
    condition_name <- paste0(TX_LABELS[tx_idx], "_T", tp)
    cat("\n--- Condition:", condition_name, "---\n")

    cond_samples <- dge$samples[dge$samples$Treatment == treatment &
                                dge$samples$Timepoint == tp, ]
    cond_expr    <- residuals[, rownames(cond_samples), drop = FALSE]
    
    # Map RNA IDs -> dbGaP VCF IDs, filter to genotyped
    sm <- data.frame(
      rna_id       = rownames(cond_samples),
      dbgap_id_vcf = gsub("\\.", "-", cond_samples$dbGaP_ID),
      stringsAsFactors = FALSE
    )
    sm <- sm[sm$dbgap_id_vcf %in% vcf_samples, ]
 
    if (nrow(sm) == 0) { warning("No genotyped samples for ", condition_name); next }
    cat(condition_name, "- genotyped:", nrow(sm), "/", nrow(cond_samples), "\n")
    
    cond_expr_filt <- cond_expr[, sm$rna_id, drop = FALSE]
    
    cat("Inverse normal transforming...\n")
    int_expr <- t(apply(cond_expr_filt, 1, function(x) {
    qqnorm(rank(x, ties.method = "random"), plot = FALSE)$x
    }))
    rownames(int_expr) <- rownames(cond_expr_filt)
    colnames(int_expr) <- colnames(cond_expr_filt)

    # PCA on INT data - extract top 30 PCs
    cat("Running PCA...\n")
    n.pc                <- min(30, ncol(int_expr) - 1)
    pca                 <- prcomp(t(int_expr), center = TRUE, scale. = FALSE)
    pc_matrix           <- pca$x[, 1:n.pc, drop = FALSE]
    colnames(pc_matrix) <- paste0("PC", 1:n.pc)
    cat("Extracted", n.pc, "PCs\n")

    # Save PC file
    pc_cond           <- t(pc_matrix)
    colnames(pc_cond) <- sm$dbgap_id_vcf
    pc_out  <- data.frame(id = rownames(pc_cond), pc_cond, check.names = FALSE)
    pc_file <- file.path(PC_DIR, paste0("PC_", condition_name, ".txt"))
    write.table(pc_out, pc_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Saved PC:", pc_file, "\n")
 
    # Save INT BED file
    int_out           <- int_expr
    colnames(int_out) <- sm$dbgap_id_vcf
 
    # FIXED 6/1/26 with TSS 
    gene_data <- anno[match(rownames(int_out), anno$ensgene), ]

    gene_data <- gene_data %>%
    mutate(
        tss  = ifelse(strand == 1, start, end - 1),
     tss  = as.integer(tss),
     tss1 = as.integer(tss + 1)
     )

    bed_int <- data.frame(
      chr     = gene_data$chr,
      start   = gene_data$tss,
      end     = gene_data$tss1,
      gene_id = rownames(int_out),
      int_out,
      check.names = FALSE
    )
    
    bed_int  <- bed_int[order(as.integer(bed_int$chr), bed_int$start), ]
    bed_file <- file.path(OUTPUT_DIR, paste0("GxP-eQTL_", condition_name, "_int.bed.gz"))
    gz_con   <- gzfile(bed_file, "w")
    writeLines(paste0("#", paste(colnames(bed_int), collapse = "\t")), gz_con)
    write.table(bed_int, gz_con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(gz_con)
    cat("Saved INT BED:", bed_file, "| Genes:", nrow(bed_int), "Samples:", ncol(int_out), "\n")
  }
}
