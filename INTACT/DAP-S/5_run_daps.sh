#!/bin/bash
# 5_run_daps.sh
# Split gene list into chunks and submit one DAP-S job per chunk per condition.

###############################################################################
# CONFIG
###############################################################################

INDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/torus_input"
OUTDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/daps_output/open_chromatin"
CHUNKDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/gene_chunks"
SCRIPTDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/DAP-S"
LOGDIR="/rs/rs_grp_gxp/RNAseq_analysis/phthalates/INTACT/daps_logs"

CHUNK_SIZE=1000
CONDITIONS=(EtOH_T6 EtOH_T24 MBP_T6 MBP_T24)

mkdir -p "${OUTDIR}" "${LOGDIR}"

###############################################################################
# 1. Build gene list and split into chunks (done once)
###############################################################################

mkdir -p "${CHUNKDIR}" "${OUTDIR}" "${LOGDIR}"

if [ ! "$(ls -A ${CHUNKDIR}/splitGene* 2>/dev/null)" ]; then
  echo "Splitting gene list into chunks of ${CHUNK_SIZE}..."
  zcat "${INDIR}/split/MBP_T6.eQTL.txt.gz" \
    | awk 'NR>1 {print $2}' \
    | sort -u \
    > "${CHUNKDIR}/gene_list.txt"

  split -l ${CHUNK_SIZE} "${CHUNKDIR}/gene_list.txt" \
        -d -a 3 "${CHUNKDIR}/splitGene"

  echo "Total genes: $(wc -l < ${CHUNKDIR}/gene_list.txt)"
  echo "Chunks created: $(ls ${CHUNKDIR}/splitGene* | wc -l)"
else
  echo "Gene chunks already exist, skipping split."
fi

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
