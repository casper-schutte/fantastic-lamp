#!/bin/bash
# This script will be used to the create the genome graph based on the reference
# genome and the edits recorded in the CSV file.

env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'
export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2

cat genomic_edits.csv | awk -F',' '{print ">homology_arm_"$1; print $8;}' | tr -d \- > homology_arms.fa
cat genomic_edits.csv | awk -F',' '{print ">ref_homology_arm_"$1; print $7;}' | tr -d \- > ref_subpaths.fa

# Combine homology arms and reference over the range of the homology arms into one FASTA file.
cat ref_subpaths.fa homology_arms.fa > ref_and_hom_arms.fa

ref=lambda_phage.fasta

# map the homology arms against the reference
minimap2 -k 19 -w 1 -cx sr $ref ref_and_hom_arms.fa >ref_and_hom_arms.paf

cat $ref ref_and_hom_arms.fa > ref+edits.fa

# induce the variation graph
seqwish -g ref+edits.gfa -s ref+edits.fa -p ref_and_hom_arms.paf -P

# sort and "chop" the graph so nodes are <256bp long (needed for vg map)
odgi build -g ref+edits.gfa -o -  | odgi sort -i - -p sYYgs -o - | odgi chop -i - -o - -c 256 | tee ref+edits.og | odgi view -i - -g >ref+edits.og.gfa

# import the graph into xg format (efficient static graph model)
vg convert -x ref+edits.og.gfa >ref+edits.og.gfa.xg

# index the graph
vg index -p -g ref+edits.og.gfa.gcsa -t 16 ref+edits.og.gfa.xg

cat Data_names.txt | while read -r line; do
  vg map -x ref+edits.og.gfa.xg -g ref+edits.og.gfa.gcsa -t 16 -% -f "$line".fastq.gz | pv -l >"$line".gaf
  python3 compare_coverage_read_info.py --gaf-path "$line".gaf --out-path "$line".tsv --og-path "ref+edits.og"
done
