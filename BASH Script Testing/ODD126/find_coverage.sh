#!/bin/bash

# Need to perform all the initial steps in commands.sh, as well as the injection
# of the reference genome over the range of the homology arms. Static graphs
# Alternatively, we could just perform copy over the files that do not need to change for each
# experiment.
# Then, the GAF file (reads aligned to the graph) is constructed for each data set. Iterative graph creation.
# Finally, the Python script will be called to perform the coverage analysis.

# 1) Make static graphs: (from commands.sh)
# extract the homology arms as FASTA
cat ODD126_homology_arms.csv | tr , '\t' | awk '{print ">homology_arm_"$1"; "$2; }' > ODD126_homology_arms.fa

# The names for the homology arms are in DesignLibraryDetails_ODD126.csv, column 1. The corresponding reference seq is
# in column BA, the 53rd column.
cat DesignLibraryDetails_ODD126.csv | awk -F',' '{print ">ref_homology_arm_"$1; print $53;}' | tr -d \- > ref_subpaths.fa

cat ref_subpaths.fa ODD126_homology_arms.fa > ODD126_ref_and_hom_arms.fa

# map the homology arms against the reference
minimap2 -k 19 -w 1 -cx sr GCA_010356925.1_ASM1035692v1_genomic.fna ODD126_ref_and_hom_arms.fa >ODD126_ref_and_hom_arms.paf

# combine the inputs to seqwish in a single file
cat GCA_010356925.1_ASM1035692v1_genomic.fna ODD126_ref_and_hom_arms.fa >yeast+edits.fa

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

# 2) Make GAF files for each data set:
# example:
# vg map -x yeast+edits.og.gfa.xg -g yeast+edits.og.gfa.gcsa -t 16 -% -f GE00001631-DOT_A07_S103_R1_001 | pv -l >GE00001631-DOT_A07_S103_R1_001.subset.gaf


# 3) Run coverage analysis with Python script:


