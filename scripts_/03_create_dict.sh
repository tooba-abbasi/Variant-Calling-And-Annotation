#!/bin/bash

# Script to create sequence dictionary using GATK

REF_FASTA="/home/proteindynamics2/NGS/ref/Homo_sapiens_assembly38.fasta"
DICT_OUTPUT="/home/proteindynamics2/NGS/ref/Homo_sapiens_assembly38.dict"

# Create sequence dictionary
gatk CreateSequenceDictionary \
    -R "$REF_FASTA" \
    -O "$DICT_OUTPUT"

echo "✅ Dictionary created at: $DICT_OUTPUT"
