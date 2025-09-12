#!/bin/bash
set -euo pipefail

# Directory where your full known‐sites VCFs live
KNOWN_DIR="$HOME/NGS/ref"
# Output directory for chr17‐only subset
OUT_DIR="$KNOWN_DIR/chr17_only"

mkdir -p "$OUT_DIR"

# Automatically find all .vcf.gz files in KNOWN_DIR
for full in "$KNOWN_DIR"/*.vcf.gz; do
  base=$(basename "$full" .vcf.gz)
  echo " › Subsetting $base → ${base}.chr17.vcf.gz"
  bcftools view -r chr17 \
      "$full" \
    -Oz -o "$OUT_DIR/${base}.chr17.vcf.gz"
  tabix -p vcf "$OUT_DIR/${base}.chr17.vcf.gz"
done

echo "✅ All known‐sites files subset to chr17 in: $OUT_DIR"
