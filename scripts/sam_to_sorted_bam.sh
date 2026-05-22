#!/usr/bin/env zsh

input_sam=$1

if [ -z "$input_sam" ]; then
    echo "Usage: $0 <input.sam>"
    exit 1
fi

samtools view -bS "$input_sam" | samtools sort -o "${input_sam:r}.sorted.bam"
samtools index "${input_sam:r}.sorted.bam"
# rm "$input_sam"

