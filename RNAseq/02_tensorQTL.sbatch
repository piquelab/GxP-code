#!/usr/bin/env bash
#SBATCH -J tQTL_cis
#SBATCH -q primary
#SBATCH -n 4
#SBATCH --mem=24G
#SBATCH -t 08:00:00
#SBATCH -o %x.%A_%a.out
#SBATCH -e %x.%A_%a.err

set -euo pipefail

################################################################################
# 02_tensorQTL.sbatch
# Purpose: Run tensorQTL cis-eQTL analysis for one condition
# Date: 02/05/2026
################################################################################

## SET USER PATHS
BASE="/rs/rs_grp_gxp/RNAseq_analysis/GxP_01092026"
PLINK_PREFIX="${BASE}/tensor/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"

# These are passed via --export from 01_tensorQTL.sh:
# COND_TSV - path to condition table for this SV iteration
# ITER_DIR - output directory for this SV iteration
# SV_LABEL - label for this iteration (e.g., SV0, SV1-5, etc.)

## LOAD ENVIRONMENT
module purge
module load plink/2.0
source activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH

python3 -m tensorqtl --help >/dev/null 2>&1 || { echo "tensorQTL not on PATH"; exit 1; }

## ARRAYED TASK
line=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==1{next} NR==n+1{print; exit}' "${COND_TSV}") || true

if [[ -z "${line:-}" ]]; then
  echo "No line for task ${SLURM_ARRAY_TASK_ID} in ${COND_TSV}" >&2
  exit 0
fi

IFS=$'\t' read -r CONDITION PHENO_BED COV_FILE <<< "${line}"

## OUTPUTS
OUT_PREFIX="${CONDITION}_cis"

echo "[$(date)] Starting tensorQTL cis for ${SV_LABEL}"
echo "SV iteration: ${SV_LABEL}"
echo "Condition   : ${CONDITION}"
echo "Genotypes   : ${PLINK_PREFIX}.[pgen|pvar|psam]"
echo "Phenotypes  : ${PHENO_BED}"
echo "Covariates  : ${COV_FILE}"
echo "Out dir     : ${ITER_DIR}"
echo "Out prefix  : ${OUT_PREFIX}"

## RUN
python3 -m tensorqtl \
  "${PLINK_PREFIX}" \
  "${PHENO_BED}" \
  "${OUT_PREFIX}" \
  --covariates "${COV_FILE}" \
  --mode cis \
  --fdr 0.1 \
  --window 100000 \
  -o "${ITER_DIR}"

echo "[$(date)] Done."
