#!/usr/bin/zsh

samtools view -h $1 \
    | awk '($6 ~ /N/) || ($1 ~ /^@/)' \
    | samtools view -bS - \
    | bamToBed -bed12
