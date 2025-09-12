#!/bin/bash

# Set variables
REF_DIR="/home/proteindynamics2/NGS/ref"
REF_PREFIX="$REF_DIR/Homo_sapiens_assembly38"
READS_DIR="/home/proteindynamics2/NGS/samples"
OUT_DIR="/home/proteindynamics2/NGS/alignments"

mkdir -p "$OUT_DIR"

# Input FASTQ files
READ1="$READS_DIR/Example_1.fastq"
READ2="$READS_DIR/Example_2.fastq"
SAMPLE="Example"

# Output BAM
SAM_OUT="$OUT_DIR/${SAMPLE}.sam"
BAM_OUT="$OUT_DIR/${SAMPLE}.bam"
SORTED_BAM="$OUT_DIR/${SAMPLE}_sorted.bam"

# Run BWA MEM
bwa mem -t 4 "$REF_PREFIX" "$READ1" "$READ2" > "$SAM_OUT"

# Convert SAM to BAM, then sort
samtools view -@ 4 -Sb "$SAM_OUT" | samtools sort -@ 4 -o "$SORTED_BAM"

# Clean up (optional)
rm "$SAM_OUT"

# Index the sorted BAM
samtools index "$SORTED_BAM"


