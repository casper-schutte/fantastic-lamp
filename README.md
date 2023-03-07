# Fantastic-lamp:

## Project Description:
This project aims to evaluate the success of genome editing by measuring the coverage of reads mapping to edited 
regions compared to the corresponding reference sequence. This is accomplished by aligning reads from the edited genome
to a genome graph built from the reference and intended edits. The pipeline can simultaneously calculate the coverage 
of multiple populations, making it an efficient tool for quantifying the success of novel editing methods or verifying 
multiple edits. The efficacy of the edits can be inferred from the output, which is a TSV file containing the list of 
intended edits, their homology coverage, and their reference coverage. 

## Installation and Dependencies:
This pipeline is initiated with the "find_coverage.sh" script from the command line and requires no explicit installation.
All the dependencies can be installed with Conda from the "environment.yaml" file.
This pipeline was developed and tested on 
Ubuntu 20.04 with Python 3.10, although earlier versions of Python may also be compatible. For more information about 
the system configuration that has been confirmed to run this pipeline correctly, please refer to the "Test.yml" file in
the /workflows directory. 


## Descriptions of steps and files used by the pipeline:

1) Homology arms (hom_arms) and the reference sequence for each homology arm (ref_hom_arms)
are extracted from DesignLibraryDetails_ODD126.withEditWindow.csv and combined into a single file: ODD126_ref_and_hom_arms.fa

2) minimap2 is used to map the hom_arms and ref_hom_arms to the reference sequence (ref_and_mt.fna) which includes both the reference 
 sequence and its mitochondrial (mtDNA) sequence. The alignment is saved as ODD126_ref_and_hom_arms.paf
3) ref_and_mt.fna and ODD126_ref_and_hom_arms.fa are combined with ODD126_augmented_CB39.fasta (this is the plasmid sequence)
to make yeast+edits.fa
4) Using yeast+edits.fa and the alignment from step 2, seqwish is used to create the variation graph (yeast+edits.gfa).
5) The graph is sorted and chopped using odgi and then converted into xg format (yeast+edits.og.gfa.xg), before finally being 
indexed -> yeast+edits.og.gfa.gcsa
6) The file Data_names.txt contains the names of the files which contain the sequencing reads. They 
are in .fastq.gz format. The script can handle paired-end reads. This can be changed in the bash script (the files need to be named appropriately)
This file is iterated over and for each line (file) the following steps (7 & 8) are executed:

7) The reads from both files are mapped onto the graph (yeast+edits.og.gfa.xg), creating "filename".gaf
8) The python script "compare_coverage.py" is called with the .gaf file and the yeast+edits.og file as input.
From the yeast+edits.og file, a dictionary is created mapping node_id to path names. From this dictionary, homology arm
and reference homology arm paths are created and edges are created from the nodes. Edges that are shared
between the ref_hom_arms and hom_arms are discarded. Then, the number of reads mapping to edges within 
hom_arms and ref_hom_arms are counted (from the .gaf file) and put into a dictionary mapping edges to the read count for that
edge. Coverage for a path calculated as the sum of the number of reads mapping to an edge in the path divided by 
the number of edges in the path. These coverages are written to a .tsv file.


