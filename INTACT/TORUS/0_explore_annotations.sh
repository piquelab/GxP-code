#!/bin/bash
# 0_explore_annotations.sh
# Inspect CENTIPEDE footprint annotation files before building TORUS inputs.

ANNOT_DIR="/rs/rs_grp_gxp/hp8265_ATAC_output/Footprinting/CENTIPEDE/Merged_AllBatch_All_Strict_Motifs_SNP_filtered"

echo "=== Total motif files ==="
ls "${ANNOT_DIR}"/*.bed.gz | wc -l

echo ""
echo "=== First 20 motif file names ==="
ls "${ANNOT_DIR}"/*.bed.gz | xargs -n1 basename | head -20

echo ""
echo "=== Column count (should be consistent across files) ==="
zcat "${ANNOT_DIR}/MA0655.1.bed.gz" | awk '{print NF}' | sort | uniq -c

echo ""
echo "=== First 3 rows of MA0655.1 (confirm column 9 = SNP ID) ==="
zcat "${ANNOT_DIR}/MA0655.1.bed.gz" | head -3 | cut -f1-10 | column -t

echo ""
echo "=== Chromosome naming (col 1) ==="
zcat "${ANNOT_DIR}/MA0655.1.bed.gz" | awk '{print $1}' | sort -u | head -10

echo ""
echo "=== Unique SNPs per motif (all files) ==="
for f in "${ANNOT_DIR}"/*.bed.gz; do
    motif=$(basename "$f" .bed.gz)
    n=$(zcat "$f" | awk '{print $9}' | sort -u | wc -l)
    echo -e "${motif}\t${n}"
done

echo ""
echo "=== Total unique SNPs across all motifs ==="
for f in "${ANNOT_DIR}"/*.bed.gz; do
    zcat "$f" | awk '{print $9}'
done | sort -u | wc -l

echo ""
echo "=== SNP ID format check (col 9, first 10 unique values from MA0655.1) ==="
zcat "${ANNOT_DIR}/MA0655.1.bed.gz" | awk '{print $9}' | sort -u | head -10

echo ""
echo "=== Check for any header lines (lines not starting with a digit or 'chr') ==="
zcat "${ANNOT_DIR}/MA0655.1.bed.gz" | awk 'NR<=5 && $1 !~ /^[0-9]/ && $1 !~ /^chr/' | head -5
echo "(empty = no header)"
