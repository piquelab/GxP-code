#!/bin/bash
# 3_run_torus.sh
# Submit one TORUS job per condition.
# Runs enrichment estimation (-est) and dumps per-gene priors (-dump_prior)
# for use in downstream DAP-G fine-mapping.

###############################################################################
# CONFIG
###############################################################################

TORUS=/wsu/el7/groups/piquelab/misc/bin/torus.static

INDIR=/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_input
OUTDIR=/rs/rs_grp_gxp/RNAseq_analysis/RNA_021826/INTACT/torus_output

# Annotation file — switch between open chromatin and CENTIPEDE by changing this line
ANNOT=${INDIR}/zzz_torus_oc_nopeaking.annot.gz   # open chromatin
# ANNOT=${INDIR}/zzz_torus.annot.gz    # CENTIPEDE

CONDITIONS=(BPA_T6 BPA_T24 EtOH_T6 EtOH_T24 H2O_T6 H2O_T24 MBP_T6 MBP_T24)

###############################################################################

mkdir -p ${OUTDIR}
 
for COND in "${CONDITIONS[@]}"; do
 
  PRIOR_DIR=${OUTDIR}/${COND}_dump_prior
 
  # Job 1: enrichment estimation
#  sbatch \
#    -q express \
#    --mem=60G \
#    --time=2-00:00:00 \
#    -N 1-1 -n 1 \
#    --job-name=torus_est_${COND} \
#    --output=${OUTDIR}/slurm_torus_est_${COND}.out \
#    --wrap "${TORUS} \
#      -d ${INDIR}/split/${COND}.eQTL.txt.gz \
#      -smap ${INDIR}/zzz_snp.map.gz \
#      -gmap ${INDIR}/zzz_gene.map.gz \
#      -annot ${ANNOT} \
#      -est > ${OUTDIR}/${COND}.est"
# 
#  echo "Submitted est: ${COND}"
#  sleep 0.5
 
  # Job 2: dump priors for DAP-G
  sbatch \
    -q express \
    --mem=60G \
    --time=2-00:00:00 \
    -N 1-1 -n 1 \
    --job-name=torus_prior_${COND} \
    --output=${OUTDIR}/slurm_torus_prior_${COND}.out \
    --wrap "${TORUS} \
      -d ${INDIR}/split/${COND}.eQTL.txt.gz \
      -smap ${INDIR}/zzz_snp.map.gz \
      -gmap ${INDIR}/zzz_gene.map.gz \
      -annot ${ANNOT} \
      -dump_prior ${PRIOR_DIR}"
 
  echo "Submitted prior: ${COND}"
  sleep 0.5
  
 done
