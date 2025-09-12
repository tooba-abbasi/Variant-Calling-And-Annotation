#!/bin/bash

# Script: Run FastQC on paired-end FASTQ files
# Project: gatk_variant_calling

SAMPLES_DIR="/home/proteindynamics2/NGS/samples"
OUTPUT_DIR="/home/proteindynamics2/NGS/fastqc_reports"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"



# Run FastQC on all FASTQ files in the samples directory
echo "🔍 Running FastQC on FASTQ files in: $SAMPLES_DIR"
fastqc "$SAMPLES_DIR"/*.fastq -o "$OUTPUT_DIR" --extract

echo "✅ FastQC completed. Reports saved in: $OUTPUT_DIR"
