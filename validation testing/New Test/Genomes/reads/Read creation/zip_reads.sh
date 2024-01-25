#!/bin/bash

# This script uses gzip to convert the fastq files into fastq.gz

for i in {1..100}; do
    gzip genome${i}_1.fastq 
    gzip genome${i}_2.fastq
done
