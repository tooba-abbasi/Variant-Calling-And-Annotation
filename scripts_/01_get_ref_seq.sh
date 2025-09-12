#!/bin/bash

# this scripts download the Human ref sequence in hg38 assembly from the gatk resouces list.

mkdir -p ref && cd ref

# Download the full GRCh38 reference with decoy sequences and no ALT contigs
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta


