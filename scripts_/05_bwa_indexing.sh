#!/bin/bash

# Define directories and reference genome path
REF="/home/proteindynamics2/NGS/ref/Homo_sapiens_assembly38.fasta"
INDEX_DIR="/home/proteindynamics2/NGS/bwa_index"

# Create directory for index if it doesn't exist
mkdir -p "$INDEX_DIR"

# Copy reference genome to the indexing directory
cp "$REF" "$INDEX_DIR"

# Change to the index directory
cd "$INDEX_DIR"

# Check if the reference fasta file exists
if [ ! -f "$REF" ]; then
    echo "Reference genome file not found!"
    exit 1
fi

# Start BWA indexing with 24 threads (adjust if needed)
echo "Starting BWA indexing with 24 threads..."

# Run BWA index command with parallelization
bwa index -p Homo_sapiens_assembly38 -a bwtsw "$REF" &

# Wait for indexing to complete
wait

echo "BWA indexing completed successfully!"
