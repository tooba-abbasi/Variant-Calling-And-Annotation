#!/bin/bash
set -euo pipefail

# Paths
REF="/home/proteindynamics2/NGS/ref/Homo_sapiens_assembly38.fasta"
WORKDIR="/home/proteindynamics2/NGS/variants"

RAW_VCF="$WORKDIR/Example.chr17.raw.vcf.gz"
SNP_VCF="$WORKDIR/Example.chr17.snps.vcf.gz"
INDEL_VCF="$WORKDIR/Example.chr17.indels.vcf.gz"

SNP_FILTERED="$WORKDIR/Example.chr17.snps.filtered.vcf.gz"
INDEL_FILTERED="$WORKDIR/Example.chr17.indels.filtered.vcf.gz"
MERGED_FILTERED="$WORKDIR/Example.chr17.filtered.vcf.gz"

mkdir -p "$WORKDIR"

echo "1) Splitting raw VCF into SNPs and INDELs..."
gatk SelectVariants -R "$REF" -V "$RAW_VCF" --select-type SNP   -O "$SNP_VCF"
gatk SelectVariants -R "$REF" -V "$RAW_VCF" --select-type INDEL -O "$INDEL_VCF"

echo "2) Hard-filtering SNPs..."
gatk VariantFiltration \
  -R "$REF" -V "$SNP_VCF" \
  --filter-name "QD_lt2"               --filter-expression "QD < 2.0" \
  --filter-name "FS_gt60"              --filter-expression "FS > 60.0" \
  --filter-name "SOR_gt3"              --filter-expression "SOR > 3.0" \
  --filter-name "MQ_lt40"              --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum_lt-12.5"    --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum_lt-8"  --filter-expression "ReadPosRankSum < -8.0" \
  -O "$SNP_FILTERED"

echo "3) Hard-filtering INDELs..."
gatk VariantFiltration \
  -R "$REF" -V "$INDEL_VCF" \
  --filter-name "QD_lt2"               --filter-expression "QD < 2.0" \
  --filter-name "FS_gt200"             --filter-expression "FS > 200.0" \
  --filter-name "ReadPosRankSum_lt-20" --filter-expression "ReadPosRankSum < -20.0" \
  -O "$INDEL_FILTERED"

echo "4) Merging filtered SNPs and INDELs..."
gatk MergeVcfs -I "$SNP_FILTERED" -I "$INDEL_FILTERED" -O "$MERGED_FILTERED"

echo "✅ Filtering complete. Final VCF: $MERGED_FILTERED"
