#!/bin/bash
# Convert TensorQTL parquet files to tab-delimited text files

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR="/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor"
INPUT_DIR="${BASE_DIR}/tensorqtl_output_cis-nominal_SV15_100kb"
OUTPUT_DIR="${BASE_DIR}/txt_files"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "[$(date)] Starting parquet to txt conversion..."
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"

# Find all parquet files
PARQUET_FILES=(${INPUT_DIR}/*.parquet)

echo "Found ${#PARQUET_FILES[@]} parquet files to convert"

# =============================================================================
# CONVERT EACH FILE
# =============================================================================
converted=0
failed=0

for parquet_file in "${PARQUET_FILES[@]}"; do
    # Get base filename without extension
    base_name=$(basename "${parquet_file}" .parquet)
    txt_file="${OUTPUT_DIR}/${base_name}.txt"
    
    echo ""
    echo "Converting: $(basename "${parquet_file}")"
    
    # Use Python to convert parquet to txt
    python3 << EOF
import pandas as pd
import sys
import os

try:
    # Read parquet file
    df = pd.read_parquet('${parquet_file}')
    
    print(f"  Shape: {df.shape[0]:,} rows x {df.shape[1]} columns")
    print(f"  Columns: {', '.join(df.columns)}")
    
    # Save as tab-delimited text file
    df.to_csv('${txt_file}', sep='\t', index=False)
    
    print(f"  Saved: $(basename "${txt_file}")")
    print(f"  Size: $(du -h "${txt_file}" | cut -f1)")
    
except Exception as e:
    print(f"  ERROR: {e}")
    sys.exit(1)
EOF

    if [[ $? -eq 0 ]]; then
        ((converted++))
        echo "  ✓ SUCCESS"
    else
        ((failed++))
        echo "  ✗ FAILED"
    fi
done

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "=============================================================="
echo "[$(date)] Conversion complete!"
echo "=============================================================="
echo "Successfully converted: ${converted} files"
echo "Failed: ${failed} files"
echo ""
echo "Output directory: ${OUTPUT_DIR}"

if [[ ${converted} -gt 0 ]]; then
    echo ""
    echo "Sample of converted files:"
    ls -lh "${OUTPUT_DIR}" | head -10
    
    echo ""
    echo "Quick preview of first converted file:"
    first_txt=$(ls "${OUTPUT_DIR}"/*.txt | head -1)
    echo "File: $(basename "${first_txt}")"
    head -5 "${first_txt}"
fi

echo ""
echo "Conversion complete!"
