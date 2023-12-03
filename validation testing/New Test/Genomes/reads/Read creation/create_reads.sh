#!/bin/bash

# This script will be used to create the reads that will be used for validation testing.

#python readSimulator.py --input genome0.fasta --simulator wgsim --simulator_path wgsim --outdir shredded_reads --iterations 1 --readlen 200 --depth 20 --opts '-e 0 -r 0 -R 0 -X 0 -h -S 5'

#wgsim -e 0.020 -r 0 -R 0 -X 0 -h -S 0 -d 200 -N 5000 -1 200 -2 200 -S 5 genome1.fasta genome1_1.fastq genome1_2.fastq

for i in {1..99}; do
    wgsim -e 0.00 -r 0 -R 0 -X 0 -h -S 0 -d 200 -N 1000 -1 200 -2 200 -S 5 genome${i}.fasta genome${i}_1.fastq genome${i}_2.fastq
done

