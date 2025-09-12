#!/bin/bash
set -euo pipefail

# Paths
REF_FASTA="$HOME/NGS/ref/Homo_sapiens_assembly38.fasta"
INPUT_BAM="$HOME/NGS/alignments/Example_sorted_with_rg_marked.bam"
OUT_DIR="$HOME/NGS/bqsr"
KNOWN_DIR="$HOME/NGS/ref/chr17_only"

DBSNP="$KNOWN_DIR/1000G_phase1.snps.high_confidence.hg38.chr17.vcf.gz"
MILLS="$KNOWN_DIR/Mills_and_1000G_gold_standard.indels.hg38.chr17.vcf.gz"

RECAL_TABLE="$OUT_DIR/recal_data.table"
POST_RECAL_TABLE="$OUT_DIR/post_recal_data.table"
RECAL_BAM="$OUT_DIR/Example_recalibrated.bam"
DIFF_PLOT="$OUT_DIR/bqsr_difference.pdf"

mkdir -p "$OUT_DIR"

# 1) First pass: build recalibration table
gatk BaseRecalibrator \
  -R "$REF_FASTA" \
  -I "$INPUT_BAM" \
  --known-sites "$DBSNP" \
  --known-sites "$MILLS" \
  -O "$RECAL_TABLE"

# 2) Apply recalibration to BAM
gatk ApplyBQSR \
  -R "$REF_FASTA" \
  -I "$INPUT_BAM" \
  --bqsr-recal-file "$RECAL_TABLE" \
  -O "$RECAL_BAM"

# 3) Second pass: build post-recalibration table on the recalibrated BAM
gatk BaseRecalibrator \
  -R "$REF_FASTA" \
  -I "$RECAL_BAM" \
  --known-sites "$DBSNP" \
  --known-sites "$MILLS" \
  -O "$POST_RECAL_TABLE"

# 4) (Optional) Generate before/after plots
gatk AnalyzeCovariates \
  -before "$RECAL_TABLE" \
  -after "$POST_RECAL_TABLE" \
  -plots "$DIFF_PLOT"

echo "✅ BQSR complete. Recalibrated BAM: $RECAL_BAM"
