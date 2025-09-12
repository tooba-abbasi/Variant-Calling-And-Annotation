#!/bin/bash
set -euo pipefail

# Paths
REF="/home/proteindynamics2/NGS/ref/Homo_sapiens_assembly38.fasta"
GVCF="/home/proteindynamics2/NGS/variants/Example.chr17.g.vcf.gz"
OUT_DIR="/home/proteindynamics2/NGS/variants"
RAW_VCF="$OUT_DIR/Example.chr17.raw.vcf.gz"

mkdir -p "$OUT_DIR"

# Joint genotyping (single-sample)
gatk GenotypeGVCFs \
  -R "$REF" \
  -V "$GVCF" \
  -O "$RAW_VCF"

echo "✅ Raw VCF generated: $RAW_VCF"
