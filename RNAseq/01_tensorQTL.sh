#!/bin/bash
################################################################################
# 01_tensorQTL.sh
# Purpose: Prepare covariate files and condition tables for iterative SV analysis
# Date: 02/05/2026
################################################################################

BASE="/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026"
TENSOR_DIR="${BASE}/tensor"
COV_BASE_DIR="${TENSOR_DIR}/covariates"

mkdir -p "${TENSOR_DIR}" "${COV_BASE_DIR}"

# Configuration
TIMEPOINTS=(6 24)
CONDITIONS=(BPA H2O MBP EtOH)
MAX_SVS=20
GENO_PCS="${TENSOR_DIR}/genotype_pcs_tensorqtl.txt"
SV_DIR="${BASE}/Norm_SVs"

################################################################################
# Sort and index BED files if needed
################################################################################

echo "Sorting and indexing BED files..."

for bed in "${BASE}"/GxP-eQTL_*_qnorm.bed.gz; do
  if [[ ! -f "$bed" ]]; then
    echo "No BED files found"
    exit 1
  fi
  
  fname=$(basename "${bed}")
  cond="${fname#GxP-eQTL_}"
  cond="${cond%_qnorm.bed.gz}"
  
  sorted_bed="${BASE}/GxP-eQTL_${cond}_qnorm.sorted.bed.gz"
  
  if [[ ! -s "$sorted_bed" ]]; then
    echo "Creating sorted BED for ${cond}..."
    (zcat "$bed" | head -1; zcat "$bed" | tail -n +2 | sort -k1,1 -k2,2n) | bgzip > "$sorted_bed"
    tabix -p bed "$sorted_bed"
  fi
done

################################################################################
# Create covariate files and condition tables for each SV iteration
################################################################################

echo "Creating covariate files and condition tables..."

for n_sv in $(seq 0 $MAX_SVS); do
  if [[ $n_sv -eq 0 ]]; then
    sv_label="SV0"
  else
    sv_label="SV1-${n_sv}"
  fi
  
  echo ""
  echo "=== Processing ${sv_label} ==="
  
  # Create directories for this iteration
  ITER_DIR="${TENSOR_DIR}/${sv_label}"
  COV_DIR="${COV_BASE_DIR}/${sv_label}"
  mkdir -p "${ITER_DIR}" "${COV_DIR}"
  
  # Create condition table for this iteration
  COND_TSV="${ITER_DIR}/tensorqtl_conditions.tsv"
  echo -e "condition\tbed\tcovariates" > "${COND_TSV}"
  
  # Process each condition-timepoint combination
  for tp in "${TIMEPOINTS[@]}"; do
    for cond in "${CONDITIONS[@]}"; do
      
      echo "  Processing ${cond}_T${tp}..."
      
      # Input files
      BED_FILE="${BASE}/GxP-eQTL_${cond}_T${tp}_qnorm.sorted.bed.gz"
      SV_FILE="${SV_DIR}/SV_${cond}_T${tp}.txt"
      
      if [[ ! -f "$BED_FILE" ]]; then
        echo "    WARNING: Missing BED file: $BED_FILE"
        continue
      fi
      
      if [[ ! -f "$SV_FILE" ]]; then
        echo "    WARNING: Missing SV file: $SV_FILE"
        continue
      fi
      
      # Output covariate file
      COV_FILE="${COV_DIR}/covariates_${cond}_T${tp}.txt"
      
      # Create covariate file using R
      Rscript - <<EOF
# Read genotype PCs
geno_pcs <- read.table("${GENO_PCS}", header = TRUE, sep = "\t", 
                       row.names = 1, check.names = FALSE)

# Read SVs for this condition
svs <- read.table("${SV_FILE}", header = TRUE, sep = "\t", 
                  row.names = 1, check.names = FALSE)

# Get samples from BED file
bed_samples <- system("zcat ${BED_FILE} | head -1 | cut -f5-", intern = TRUE)
bed_samples <- unlist(strsplit(bed_samples, "\t"))

# Filter genotype PCs to BED samples (keep top 3 PCs)
geno_pcs_filt <- geno_pcs[1:3, bed_samples, drop = FALSE]

if ($n_sv == 0) {
  # No SVs, just genotype PCs
  cov_combined <- geno_pcs_filt
} else {
  # Subset SVs to first n_sv and to BED samples
  svs_filt <- svs[1:${n_sv}, bed_samples, drop = FALSE]
  
  # Combine genotype PCs and SVs
  cov_combined <- rbind(geno_pcs_filt, svs_filt)
}

# Add id column
cov_out <- data.frame(id = rownames(cov_combined), cov_combined, check.names = FALSE)

# Write covariate file
write.table(cov_out, "${COV_FILE}", sep = "\t", quote = FALSE, row.names = FALSE)

cat("    Covariates:", nrow(cov_out), "x", ncol(cov_out) - 1, "\n")
EOF
      
      # Add to condition table
      echo -e "${cond}_T${tp}\t${BED_FILE}\t${COV_FILE}" >> "${COND_TSV}"
      
    done
  done
  
  # Count conditions for this iteration
  N=$(( $(wc -l < "${COND_TSV}") - 1 ))
  echo "  Condition table: ${COND_TSV} (${N} conditions)"
  
  # Submit job array for this iteration
  if [[ $N -gt 0 ]]; then
    echo "  Submitting job array for ${sv_label}..."
    sbatch --array=1-"$N" --export=ALL,COND_TSV="${COND_TSV}",ITER_DIR="${ITER_DIR}",SV_LABEL="${sv_label}" 02_tensorQTL.sbatch
  else
    echo "  WARNING: No conditions to process for ${sv_label}"
  fi
  
done

echo ""
echo "=== Complete ==="
echo "Generated covariate files and submitted jobs for SV iterations 0-${MAX_SVS}"
