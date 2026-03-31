#!/bin/bash
# 9_run_enloc.sh
# Step 1: Run TORUS on GWAS z-scores to get per-SNP GWAS PIPs
# Step 2: Run fastENLOC per condition x trait to get gene-level GLCPs

###############################################################################
# CONFIG
###############################################################################

TORUS=/wsu/el7/groups/piquelab/misc/bin/torus.static
FASTENLOC=/wsu/home/groups/piquelab/fastenloc.static

GWASDIR=/nfs/rprdata/julong/sc_multiome/gwas/gwas_prepare/gwas_imputefile
INDIR=/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/enloc_input
OUTDIR=/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/enloc_output
LOGDIR=/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/enloc_logs

TOTAL_VARIANTS=2303524

TRAITS=(CARDIoGRAM_C4D_CAD)
CONDITIONS=(BPA_T6 BPA_T24 MBP_T6 MBP_T24)

###############################################################################

mkdir -p ${INDIR} ${LOGDIR}

###############################################################################
# Step 1: TORUS on GWAS (once per trait)
###############################################################################

for TRAIT in "${TRAITS[@]}"; do

  PIPFILE="${INDIR}/${TRAIT}.pip.gz"
  if [ -f "${PIPFILE}" ]; then
    echo "GWAS PIP already exists for ${TRAIT}, skipping TORUS step."
    continue
  fi

  GWAS_TORUS="${INDIR}/${TRAIT}_torus.txt.gz"
  if [ ! -f "${GWAS_TORUS}" ]; then
    sbatch \
      -q reqp \
      --mem=20G \
      --time=2:00:00 \
      -N 1-1 -n 1 \
      --job-name=gwas_prep_${TRAIT} \
      --output=${LOGDIR}/slurm_gwas_prep_${TRAIT}.out \
      --wrap="zcat ${GWASDIR}/${TRAIT}_impute_gwas.txt.gz | \
        awk 'NR>1 {print \$1, \"NA\", \$4}' | \
        gzip > ${GWAS_TORUS}"
    echo "Submitted GWAS prep: ${TRAIT}"
  fi

  sbatch \
    -q reqp \
    --mem=30G \
    --time=2-12:00:00 \
    -N 1-1 -n 1 \
    --job-name=torus_gwas_${TRAIT} \
    --output=${LOGDIR}/slurm_torus_gwas_${TRAIT}.out \
    --wrap="${TORUS} --load_zval \
      -d ${GWAS_TORUS} \
      -dump_pip ${INDIR}/${TRAIT}.pip && \
      gzip ${INDIR}/${TRAIT}.pip"

  echo "Submitted TORUS GWAS PIP: ${TRAIT}"
  sleep 0.5

done

###############################################################################
# Step 2: fastENLOC — one job per condition x trait
###############################################################################

for TRAIT in "${TRAITS[@]}"; do
  for COND in "${CONDITIONS[@]}"; do

    mkdir -p ${OUTDIR}/${TRAIT}

    EQTL_VCF="${INDIR}/${COND}_enloc.eqtl.vcf.gz"
    GWAS_PIP="${INDIR}/${TRAIT}.pip.gz"
    PREFIX="${OUTDIR}/${TRAIT}/${COND}"

    if [ -f "${PREFIX}.enloc.gene.out" ]; then
      echo "Skipping ${TRAIT} x ${COND} (already done)"
      continue
    fi

    sbatch \
      -q primary \
      --mem=60G \
      --time=12:00:00 \
      -N 1-1 -n 4 \
      --job-name=enloc_${TRAIT}_${COND} \
      --output=${LOGDIR}/slurm_enloc_${TRAIT}_${COND}.out \
      --wrap="${FASTENLOC} \
        -eqtl ${EQTL_VCF} \
        -gwas ${GWAS_PIP} \
        -total_variants ${TOTAL_VARIANTS} \
        -thread 4 \
        -prefix ${PREFIX}"

    echo "Submitted fastENLOC: ${TRAIT} x ${COND}"
    sleep 0.5

  done
done

echo "All jobs submitted."
