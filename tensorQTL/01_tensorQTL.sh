#!/bin/bash

BASE="/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping"
mkdir -p "${BASE}/tensor"

## RUN VCF TO PLINK
module load plink/2.0

VCF="${BASE}/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz"
OUT_PREFIX="${BASE}/tensor/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"

# Make PLINK 2.0 files: .pgen/.pvar/.psam (autosomes only; force biallelic IDs if missing)
plink2 \
  --vcf "${VCF}" \
  --max-alleles 2 \
  --set-missing-var-ids @:#:\$r:\$a \
  --chr 1-22 \
  --make-pgen \
  --out "${OUT_PREFIX}"

echo "Done. PLINK2 prefix: ${OUT_PREFIX}"

## BUILD CONDITION TABLE (updated for tensorQTL)

OUT_TSV="${BASE}/tensor/tensorqtl_conditions.tsv"

shopt -s nullglob
beds=( "${BASE}"/GxP-eQTL_*_qnorm_genotypePC_COMBATNone_corrected.bed.gz )

if (( ${#beds[@]} == 0 )); then
  echo "No phenotype BEDs found. Check the directory/patterns." >&2
  exit 1
fi

echo -e "condition\tbed\tcovariates" > "${OUT_TSV}"

for bed in "${beds[@]}"; do
  fname=$(basename "${bed}")
  # Extract <COND> between 'GxP-eQTL_' and '_qnorm_genotypePC_COMBATNone_corrected.bed.gz'
  cond="${fname#GxP-eQTL_}"
  cond="${cond%_qnorm_genotypePC_COMBATNone_corrected.bed.gz}"

  # sorted BED filename
  sorted_bed="${BASE}/GxP-eQTL_${cond}_qnorm_genotypePC_COMBATNone_corrected.sorted.bed.gz"

# Check if sorted BED exists; if not, create it
 if [[ ! -s "$sorted_bed" ]]; then
    echo "Creating sorted BED with header for ${cond}..."

    # Sort the BED file while preserving the original header
    (zcat "$bed" | head -1; zcat "$bed" | tail -n +2 | sort -k1,1 -k2,2n) | gzip > "$sorted_bed"
  fi

  cov="${BASE}/covariates/GxP-eQTL_${cond}-SV1-15.covariates-COMBATNone-FastQTL.txt"
  if [[ ! -s "${cov}" ]]; then
    echo "WARNING: Missing covariates for ${cond}: ${cov}" >&2
    continue
  fi

  echo -e "${cond}\t${sorted_bed}\t${cov}" >> "${OUT_TSV}"
done

echo "Wrote ${OUT_TSV}"
wc -l "${OUT_TSV}"

# Count data lines (minus header)
N=$(( $(wc -l < "${OUT_TSV}") - 1 ))
echo "N = $N"


## SUBMIT WITH ARRAY SIZED TO N
 sbatch --array=1-"$N" 02_tensorqtl.sbatch
