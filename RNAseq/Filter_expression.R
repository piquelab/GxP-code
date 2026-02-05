#!/usr/bin/env Rscript
################################################################################
# Filter_expression.R
# Purpose: Separate by time point, filter for expressed genes, TMM normalize, filter for autosomal coding genes, voom normalize
# Date: 02/05/2026
################################################################################
library(edgeR)
library(preprocessCore)
library(annotables)
library(dplyr)
library(limma)

# Configure inputs
base <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026"
input <- file.path(base, "GxP_SamplesRemoved_01182026.RData")
output <- file.path(base, "GxP_Filtered_02052026.RData")
timepoints <- c(6, 24)

# Load Data
load(input)

# Gene annotations
anno <- grch38 %>% 
  dplyr::select(ensgene, chr, start, end, biotype) %>% 
  dplyr::filter(biotype == "protein_coding", chr %in% c(1:22)) %>% 
  unique()
anno <- anno[!duplicated(anno$ensgene), ]

# Storage for results
filtered_data <- list()

# Process each timepoint
for (tp in timepoints) {
  cat("Processing timepoint:", tp, "\n")
  
  colData <- c2[c2$Timepoint == tp, ]
  counts <- m2[, rownames(colData)]
  
  ### Clean gene IDs
  gene_ids_clean <- sub("\\..*", "", rownames(counts))
  rownames(counts) <- gene_ids_clean
  
  ### Count matrix for timepoint samples and sample metadata
  dge <- DGEList(counts = counts, samples = colData)
  
  ### Apply edgeR expression filtering
  filt <- filterByExpr(dge, design = model.matrix(~Treatment, data = colData))
  dge <- dge[filt, ]
  
  ### TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")
  
  ### Filter to protein-coding autosomal
  keep <- rownames(dge) %in% anno$ensgene
  dge <- dge[keep, ]
  cat("Genes after filtering:", nrow(dge), "\n")
  
  ### Set treatment levels
  dge$samples$Batch_tx <- factor(dge$samples$Batch_tx)
  dge$samples$Sequencer <- factor(dge$samples$Sequencer)
  dge$samples$Treatment <- factor(dge$samples$Treatment, 
                                   levels = c("H2O", "EtOH", "BPA_100nM", "MBP_500nM"))
  dge$samples$dbGaP_ID <- make.names(dge$samples$dbGaP_ID)
  
  ### Voom transformation
  f <- as.formula("~Treatment + dbGaP_ID + trimmed_dClean.dFastq")
  design <- model.matrix(f, data = dge$samples)
  print(f)
  print("starting voom")
  transformed_data <- voom(dge, design = design, plot = FALSE)
  print("voom finished")
  
  ### Store results
  filtered_data[[paste0("tp", tp)]] <- list(
    dge = dge,
    voom = transformed_data,
    design = design
  )
}

# Save output
save(filtered_data, anno, file = output)
cat("Saved to:", output, "\n")
