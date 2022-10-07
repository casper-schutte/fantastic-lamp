# Fantastic-lamp:

Project description:
This project aims to evaluate the success of making edits to a genome by measuring the coverage of reads mapping to 
 edited regions compared to the reference sequence for that same region. This is achieved by aligning reads from the 
 edited genome to a genome graph constructed from the reference and the intended edits. 

Steps and file descriptions:

1) Homology arms and the reference sequence corresponding to the range of each homology arm 
are extracted from DesignLibraryDetails_ODD126.withEditWindow.csv. They are then concatenated into a single file: ODD126_ref_and_hom_arms.fa

2) minimap2 is used to map the hom_arms and ref_hom_arms to the reference: ref_and_mt.fna. This file contains the reference 
 sequence with its mitochondrial (mtDNA) sequence. Alignment saved as ODD126_ref_and_hom_arms.paf
3) ref_and_mt.fna and ODD126_ref_and_hom_arms.fa are combined with ODD126_augmented_CB39.fasta (this is the plasmid sequence)
to make yeast+edits.fa
4) Using yeast+edits.fa and the alignment from step 2, seqwish is used to create the variation graph (yeast+edits.gfa).
5) The graph is sorted and chopped using odgi and then converted into xg format (yeast+edits.og.gfa.xg), before finally being 
indexed -> yeast+edits.og.gfa.gcsa
6) The file Data_names.txt contains the names of the files which contain the sequencing reads. They 
are in .fastq.gz format. There are 2 files for each experiment, read pair 1 and 2 (R1 and R2).
This file is iterated over and for each line (file) the following steps are executed:

7) The reads from both files are mapped onto the graph (yeast+edits.og.gfa.xg), creating "filename".gaf
8) The python script "compare_coverage.py" is called with the .gaf file and the yeast+edits.og file as input.
From the yeast+edits.og file, a dictionary is created mapping node_id to path names. From this dictionary, homology arm
and reference homology arm paths are created and edges are created from the nodes. Edges that are shared
between the ref_hom_arms and hom_arms are discarded. Then, the number of reads mapping to edges within 
hom_arms and ref_hom_arms are counted (from the .gaf file) and put into a dictionary mapping edges to the read count for that
edge. Coverage for a path calculated as the sum of the number of reads mapping to an edge in the path divided by 
the number of edges in the path. These coverages are written to a .tsv file.
