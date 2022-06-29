#!/bin/bash

# The commands below are to ensure that odgi works correctly in the python script.
env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'
export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2
#export LD_PRELOAD=/home/ec2-user/.conda/pkgs/libjemalloc-5.2.1-h9c3ff4c_6/lib/libjemalloc.so.2


# This script requires as input:
# 1) a design library csv file (DesignLibraryDetails_ODD126.csv)
# 2) a reference genome with mtDNA sequence added (ref_and_mt.fna)
# 3) the python script compare_coverage.py (in the same directory as this script)
# 4) a text document (data_names.txt) containing the names of the files containing the reads, without the .fastq.gz
# extension eg: "GE00001631-DOT_H11_S191_R2_001.fastq.gz" should be "GE00001631-DOT_H11_S191_R2_001"
# 5) now also requires a fasta file containing the sequences of the plasmids used in the experiment.
# The output is in the same folder as the original files,

# 1) Make static graphs: (from commands.sh)
# extract the homology arms as FASTA
# The names for the homology arms are in DesignLibraryDetails_ODD126.csv, column 1. The corresponding reference seq is
# in column BA, the 53rd column, and the corresponding homology_arm edit sequence is in column 54.
# Note for future use: can if the design library is different, create a new variable containing the appropriate column
# numbers.
cat DesignLibraryDetails_ODD126.csv | awk -F',' '{print ">homology_arm_"$1; print $54;}' | tr -d \- > ODD126_homology_arms.fa
cat DesignLibraryDetails_ODD126.csv | awk -F',' '{print ">ref_homology_arm_"$1; print $53;}' | tr -d \- > ref_subpaths.fa

# Combine homology arms and reference over the range of the homology arms into one FASTA file.
cat ref_subpaths.fa ODD126_homology_arms.fa > ODD126_ref_and_hom_arms.fa

# map the homology arms against the reference
minimap2 -k 19 -w 1 -cx sr ref_and_mt.fna ODD126_ref_and_hom_arms.fa >ODD126_ref_and_hom_arms.paf

# combine the inputs to seqwish in a single file (and add plasmid sequences)
cat ref_and_mt.fna ODD126_ref_and_hom_arms.fa ODD126_augmented_CB39.fasta >yeast+edits.fa

# induce the variation graph
seqwish -g yeast+edits.gfa -s yeast+edits.fa -p ODD126_ref_and_hom_arms.paf -P

# sort and "chop" the graph so nodes are <256bp long (needed for vg map)
odgi build -g yeast+edits.gfa -o -  | odgi sort -i - -p sYYgs -o - | odgi chop -i - -o - -c 256 | tee yeast+edits.og | odgi view -i - -g >yeast+edits.og.gfa

# import the graph into xg format (efficient static graph model)
vg convert -x yeast+edits.og.gfa >yeast+edits.og.gfa.xg

# index the graph
vg index -p -g yeast+edits.og.gfa.gcsa -t 16 yeast+edits.og.gfa.xg
# The files generated so far can be used for all the data sets, as they do not change. Only the reads from the
# experiment change.

# 2) Make GAF files and run python script for each data set:

cat Data_names.txt | while read -r line; do
  vg map -x yeast+edits.og.gfa.xg -g yeast+edits.og.gfa.gcsa -t 16 -% -f "$line".fastq.gz | pv -l >"$line".gaf
  python3 compare_coverage.py --gaf-path "$line".gaf --out-path "$line" --og-path "yeast+edits.og"
done

echo "Done!"
