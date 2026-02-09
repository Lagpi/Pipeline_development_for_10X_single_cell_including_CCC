#!/bin/bash

#==============================================
# 00_starsolo.sh
# 
# STARsolo processing module
# Takes nextflow config parameters and processes single-cell data
#==============================================

set -euo pipefail

# Parameters from nextflow config
GENOME_DIR="$1"
FASTQ_R1="$2"
FASTQ_R2="$3"
SAMPLE_NAME="$4"
OUTPUT_DIR="$5"
WHITELIST="$6"
THREADS="$7"
SOLO_UMI_LEN="$8"
SOLO_CB_LEN="$9"
SOLO_UMI_START="${10}"
MEMORY_LIMIT="${11}"

# Setup logging
LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_starsolo.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "STARsolo Module: Processing sample $SAMPLE_NAME"
echo "Config parameters:"
echo "  UMI length: $SOLO_UMI_LEN"
echo "  CB length: $SOLO_CB_LEN" 
echo "  UMI start: $SOLO_UMI_START"
echo "  Threads: $THREADS"

# Check if output already exists
STARSOLO_OUTPUT="$OUTPUT_DIR/${SAMPLE_NAME}_Solo.out/Gene/filtered/barcodes.tsv"
if [[ -f "$STARSOLO_OUTPUT" ]]; then
    echo "STARsolo output already exists, skipping..."
    exit 0
fi

# Create temporary directory
TMP_DIR="$OUTPUT_DIR/${SAMPLE_NAME}_tmp"
mkdir -p "$TMP_DIR"

start_time=$(date +%s)

echo "Running STARsolo alignment..."

# Run STARsolo with config parameters
STAR \
    --runMode alignReads \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$FASTQ_R2" "$FASTQ_R1" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$TMP_DIR/${SAMPLE_NAME}_" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --soloBarcodeReadLength 0 \
    --soloFeatures Gene \
    --runThreadN "$THREADS" \
    --soloUMIlen "$SOLO_UMI_LEN" \
    --soloCBlen "$SOLO_CB_LEN" \
    --soloUMIstart "$SOLO_UMI_START" \
    --limitBAMsortRAM "$MEMORY_LIMIT" \
    --soloCellFilter EmptyDrops_CR \
    --outTmpDir "$TMP_DIR/STARtmp"

# Check success
if [[ ! -f "$TMP_DIR/${SAMPLE_NAME}_Solo.out/Gene/filtered/barcodes.tsv" ]]; then
    echo "Error: STARsolo failed"
    exit 1
fi

# Move results
mv "$TMP_DIR/${SAMPLE_NAME}_Solo.out" "$OUTPUT_DIR/"

# Compress Solo.out files
echo "Compressing Solo.out matrix files..."
SOLO_DIR="$OUTPUT_DIR/${SAMPLE_NAME}_Solo.out/Gene"

for dir in raw filtered; do
    if [[ -d "$SOLO_DIR/$dir" ]]; then
        echo "  Compressing $dir directory..."
        for file in "$SOLO_DIR/$dir"/*.mtx "$SOLO_DIR/$dir"/*.tsv; do
            if [[ -f "$file" && ! -f "${file}.gz" ]]; then
                gzip "$file"
                echo "    Compressed: $(basename $file)"
            fi
        done
    fi
done

echo "Compression completed"

# Cleanup
rm -rf "$TMP_DIR"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "STARsolo completed successfully in ${runtime}s"
echo "Output: $OUTPUT_DIR/${SAMPLE_NAME}_Solo.out/"