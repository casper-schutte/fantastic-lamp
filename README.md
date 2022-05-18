# fantastic-lamp:
(placeholder name)

This script is still under development. 

# Usage:
The main script is find_coverage.sh. It requires a design library CSV file, a reference genome (.fa or .fna), a text 
file containing the names of the data sets (WITHOUT the file extensions), and data sets of reads (fastq.gz) to be aligned. 
The script creates and indexes a graph based on the reference genome and homology arms. Then, for each data set, it aligns the reads to the graph, 
outputs a .gaf file, and then calculates the coverage of the homology arm and reference genome subpaths (which outputs a .tsv file). 

The output is a .tsv file containing hom_arm_name, hom_arm_coverage, ref_subpath_coverage, #_of_hom_arm_edges, 
count_for_hom_arm, #_of_ref_subpath_edges, and count_for_ref_subpath as columns.
