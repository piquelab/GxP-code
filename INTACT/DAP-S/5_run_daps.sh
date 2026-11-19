#!/bin/bash
# 5_run_daps.sh

###############################################################################
# CONFIG
###############################################################################

OUTDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/daps_output/tss"
CHUNKDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/gene_chunks"
SCRIPTDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/DAP-S"
LOGDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/daps_logs"

CONDITIONS=(EtOH_T6 EtOH_T24 MBP_T6 MBP_T24)

mkdir -p "${OUTDIR}" "${LOGDIR}"

###############################################################################
# Submit one job per chunk per condition
###############################################################################

for COND in "${CONDITIONS[@]}"; do
  mkdir -p "${OUTDIR}/${COND}"

  for CHUNK in "${CHUNKDIR}"/splitGene*; do
    CHUNK_ID=$(basename "${CHUNK}" | sed 's/splitGene//')

    OUTFILE="${OUTDIR}/${COND}/${COND}_chunk${CHUNK_ID}_pips.txt.gz"
    if [ -f "${OUTFILE}" ]; then
      echo "Skipping ${COND} chunk ${CHUNK_ID} (already done)"
      continue
    fi

    sbatch \
      -q primary \
      --mem=32G \
      --time=10:00:00 \
      -N 1-1 -n 1 \
      --job-name="daps_${COND}_${CHUNK_ID}" \
      --constraint="avx2" \
      --output="${LOGDIR}/slurm_daps_${COND}_${CHUNK_ID}.out" \
      --wrap="module load plink/2.0; Rscript ${SCRIPTDIR}/5_daps_worker_tss.R ${COND} ${CHUNK}"

    sleep 0.2
  done

  echo "Submitted all chunks for: ${COND}"
done

echo "Done submitting."
