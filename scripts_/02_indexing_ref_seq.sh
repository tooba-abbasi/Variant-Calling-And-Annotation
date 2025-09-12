#!/bin/bash

# ref sequence indexing

REF_SEQ='/home/proteindynamics2/NGS/ref/Homo_sapiens_assembly38.fasta'

samtools faidx ${REF_SEQ}

