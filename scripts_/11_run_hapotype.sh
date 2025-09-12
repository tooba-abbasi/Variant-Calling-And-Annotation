x#!/bin/bash
set -euo pipefail

# ─── Paths ────────────────────────────────────────────────────────────────
REF="/home/proteindynamics2/NGS/ref/Homo_sapiens_assembly38.fasta"
INPUT_BAM="/home/proteindynamics2/NGS/bqsr/Example_recalibrated.bam"
OUT_DIR="/home/proteindynamics2/NGS/variants"
OUTPUT_GVCF="$OUT_DIR/Example.chr17.g.vcf.gz"

# ─── Create output directory if needed ────────────────────────────────────
mkdir -p "$OUT_DIR"

# ─── Run HaplotypeCaller on chr17 only ────────────────────────────────────
gatk HaplotypeCaller \
  -R "$REF" \
  -I "$INPUT_BAM" \
  -O "$OUTPUT_GVCF" \
  --emit-ref-confidence GVCF \
  --standard-min-confidence-threshold-for-calling 30.0 \
  -L chr17

echo "✅ chr17‐restricted GVCF generated: $OUTPUT_GVCF"
