# fantastic-lamp:
(placeholder name)

This script is still under development. 

# Usage:
The script takes as input a graph of a reference genome and homology arms (in .og format), and a .gaf file containing the alignment of reads to the graph. 
The homology arm paths are identified by starting with "homology_arm_". To calculate the coverage of homology arm edits, the paths of the reference genome 
over the range of the homology arms must also be added to the .og graph. They must start with "ref_homology_arm_". The two input files are selected by 
renaming the variables og_path and gaf_path. The output is a .tsv file containing hom_arm_name, hom_arm_coverage, ref_subpath_coverage, #_of_hom_arm_edges, 
count_for_hom_arm, #_of_ref_subpath_edges, and count_for_ref_subpath as columns. The output filename is specified with the "name_for_output" variable.

