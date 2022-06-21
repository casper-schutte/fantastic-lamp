#!/bin/bash

set -ex

cat DesignLibraryDetails_ODD126.csv | tail -n +2 | awk -F',' '{print ">homology_arm_"$1; print $54;}' | tr -d \- > ODD126_homology_arms.fa
cat DesignLibraryDetails_ODD126.csv | tail -n +2 | awk -F',' '{print ">ref_homology_arm_"$1; print $53;}' | tr -d \- > ref_subpaths.fa

# Combine homology arms and reference over the range of the homology arms into one FASTA file.
cat ref_subpaths.fa ODD126_homology_arms.fa > ref_and_hom_arms.fa

# map the homology arms against the reference
minimap2 -k 19 -w 1 -cx sr GCA_016858175.1_ASM1685817v1_genomic.fna ref_and_hom_arms.fa > ref_and_hom_arms.paf
# map the plasmids against the reference

# combine the inputs to seqwish in a single file
cat GCA_016858175.1_ASM1685817v1_genomic.fna ODD126_augmented_CB39.fasta ref_and_hom_arms.fa > yeast+edits+plasmid.fa

# induce the variation graph
seqwish -g yeast+edits+plasmid.gfa -s yeast+edits+plasmid.fa -p ref_and_hom_arms.paf -P

# sort and "chop" the graph so nodes are <256bp long (needed for vg map)
odgi build -g yeast+edits+plasmid.gfa -o -  | odgi sort -i - -p sYYgs -o - | odgi chop -i - -o - -c 256 | tee yeast+edits+plasmid.og | odgi view -i - -g >yeast+edits+plasmid.og.gfa

# import the graph into xg format (efficient static graph model)
vg convert -x yeast+edits+plasmid.og.gfa >yeast+edits+plasmid.og.gfa.xg

# index the graph
vg index -p -g yeast+edits+plasmid.og.gfa.gcsa -t 16 yeast+edits+plasmid.og.gfa.xg
# The files generated so far can be used for all the data sets, as they do not change. Only the reads from the
# experiment change.
