#!/usr/bin/env bash
set -euo pipefail

# CONFIGURATION
ANNOVAR=~/tools/annovar                     # path to your ANNOVAR installation
DBDIR="${ANNOVAR}/humandb"                  # path to your humandb directory
INDIR=variants                              # where your VCFs live
OUTPREFIX=Example.chr17                     # prefix for all outputs
VCF="${INDIR}/${OUTPREFIX}.filtered.vcf.gz" # your final, hard‐filtered VCF

# 1) Download refGene if missing (required for functional annotation) :contentReference[oaicite:0]{index=0}
if [ ! -s "${DBDIR}/hg38_refGene.txt" ]; then
  echo ">> Fetching hg38 refGene database..."
  perl ${ANNOVAR}/annotate_variation.pl \
       -buildver hg38 \
       -downdb \
       -webfrom annovar \
       refGene \
       ${DBDIR}
fi

# 2) Convert VCF to ANNOVAR input (avinput) :contentReference[oaicite:1]{index=1}
echo ">> Converting VCF → ANNOVAR input..."
${ANNOVAR}/convert2annovar.pl \
     -format vcf4 \
     -allsample \
     -withfreq \
     -includeinfo \
     ${VCF} > ${INDIR}/${OUTPREFIX}.avinput

# 3) Run table_annovar.pl with all protocols
echo ">> Running table_annovar.pl annotations..."
${ANNOVAR}/table_annovar.pl \
     ${INDIR}/${OUTPREFIX}.avinput \
     ${DBDIR} \
     -buildver hg38 \
     -out ${INDIR}/${OUTPREFIX}.annotated \
     -remove \
     -protocol refGene,avsnp150,gnomad211_genome,clinvar_20221231,cosmic70,dbnsfp42a \
     -operation g,f,f,f,f,f \
     -nastring . \
     -polish \
     -otherinfo

echo "✅ ANNOVAR annotation complete."
echo "   Multianno output: ${INDIR}/${OUTPREFIX}.annotated.hg38_multianno.txt"
