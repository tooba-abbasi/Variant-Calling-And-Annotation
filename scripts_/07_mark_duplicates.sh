#!/bin/bash

# Simple script to mark duplicates in a BAM file

# Input BAM file and output BAM file
INPUT_BAM="$1"
OUTPUT_BAM="$2"

# Check if input BAM file is provided
if [[ -z "$INPUT_BAM" || -z "$OUTPUT_BAM" ]]; then
    echo "Usage: $0 <input_bam> <output_bam>"
    exit 1
fi

# Run GATK MarkDuplicates to mark duplicates in the BAM file
gatk MarkDuplicates \
    -I "$INPUT_BAM" \
    -O "$OUTPUT_BAM" \
    -M "${OUTPUT_BAM%.bam}.metrics.txt" \
    --CREATE_INDEX true

# Check if the command was successful
if [[ $? -eq 0 ]]; then
    echo "Duplicates marked successfully in $OUTPUT_BAM"
else
    echo "Error: Failed to mark duplicates"
    exit 1
fi
