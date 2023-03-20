---
title: 'FantasticLamp: a genome graph based pipeline for calculating the efficacy of genomic edits'
bibliography: references.bib
tags:
  - Python
  - Bioinformatics
  - Genome graphs
  - Genome editing
authors:
  - name: Casper J. H. Schutte
    orcid: 0000-0003-4245-6842
    equal-contrib: false
    affiliation: 1
    corresponding: false
  - name: Ian T. Fiddes
    orcid: 0000-0002-1580-7443
    equal-contrib: false
    affiliation: 2
    corresponding: false
  - name: Erik Garrison
    orcid: 0000-0003-3821-631X
    equal-contrib: false
    affiliation: 3
    corresponding: true
affiliations:
  - name: Stellenbosch University, South Africa
    index: 1
  - name: Inscripta Inc
    index: 2
  - name: The University of Tennessee Health Science Center, United States
    index: 3
date: 2 March 2023
---
# Summary

Accurately calculating the efficacy of genomic edits is crucial to understanding the performance of the editing techniques, in order to optimize methods and improve the success rate of the edits.
Additionally, by understanding the success rate of genomic edits, researchers can identify any potential problems or limitations of the techniques and work to overcome them.
This can help to improve the accuracy and precision of the editing methods, which is essential for many applications, such as creating genetically engineered cells for therapeutic purposes [@doudna2014new], understanding gene functions [@liu2019casx], and studying genetic diversity in a population of cells[@hsu2014development].
FantasticLamp is an open source pipeline for calculating the efficacy of genomic edits performed on multiple populations of cells.
It constructs a variation graph from the reference genome and edit template sequences that includes both edited and unedited sequences as paths, and combinations of edits as potential walks through the graph.
It then maps reads from the edited populations onto this graph in order to calculate the coverage of the edited sequences compared to the unedited sequences.
This pipeline aims to provide a quantitative measure of the success of each genomic edit without bias towards the reference or particular single editing constructs.

# Statement of need

In the field of genome engineering, researchers often perform genomic edits in cells, such as CRISPR/Cas9, TALEN, and ZNF-based systems, to study gene functions and to create new cell lines [@gaj2013zfn].
However, these edits may not always be successful, and it may be challenging to identify and quantify the success rate of these edits in large-scale data sets [@guell2014genome], [@van2020delivery].
Compared to linear alignment methods, genome graphs can provide a more accurate and comprehensive view [@garrison2018variation] of the relationships between sequences [@paten2017genome], and by using this method, one can identify and quantify the success rate of genomic edits.
When assessing the success rate of genomic edits, a linear alignment method may introduce reference bias and overlook the complexity of edits, particularly when multiple nearby edits and edit state mixtures are present.
Linear approaches that attempt to consider all possible edit states as individual sequences in order to accurately represent complexity can become overly complex themselves and prone to errors [@huang2013short], [@mun2021leviosam].
Alternatively, a genome graph approach that attempts the same thing can provide a more complete and accurate perspective of sequence relationships [@eggertsson2017graphtyper], while avoiding reference bias and capturing all allele state mixtures present within the targeted edits, even in the case of overlapped edits.
The cost of this increased accuracy and the ability to represent complexity is that using genome graphs is more computationally intensive than linear alignments methods [@rakocevic2019fast]. 

FantasticLamp is a pipeline that consists of a bash script that can be initiated from the command line.
The bash script calls a Python script written for this pipeline and the following bioinformatics tools: vg, which constructs genome graphs to represent genetic variation and facilitates efficient variant analysis [@garrison2018variation]; odgi, which optimizes the representation of sequence graphs for scalable genome analysis and visualization [@guarracino2022odgi]; minimap2, which rapidly aligns long sequencing reads to a reference sequence, enabling efficient variant calling and structural variation analysis [@li2018minimap2];
and seqwish, which converts sequences and independently generated alignments between them into a variation graph [@garrison2023unbiased] -- here used to generate the reference graph target for alignment.
FantasticLamp was designed to calculate the coverage of genomic edits in multiple populations of cells, simultaneously.
The pipeline can handle both paired-end and single-end reads.


Given a reference genome and reads sequenced from multiple edited populations, the pipeline uses a design library CSV file, which contains the intended edits and reference sequences at the intended edit sites, to construct a genome graph.
The genome graph is made up of the reference genome, the intended edit sequences ("homology arms", or "edit homology arms"), and the reference sequences at the edit sites ("reference homology arms").
The construction of the graph involves taking the homology arms and reference homology arms and mapping them to the reference genome using minimap2 [@li2018minimap2]. 
A variation graph is then induced by seqwish [@garrison2023unbiased] that represents the relationships between the reference genome and both groups of homology arms (reference and edit).
The inclusion of the reference homology arms is what allows the pipeline to calculate the relative coverage (and thus the edit efficiency) by comparing the edit homology arm coverage to the reference homology arm coverage.
Next, the graph is constructed using odgi [@guarracino2022odgi] by first "chopping" the nodes into segments smaller than 256 base pairs and then sorting the graph. These steps are necessary in order for the following steps to work.
The graph is converted into a format that is more efficient and finally indexed, both steps using vg [@garrison2018variation].
The pipeline then maps reads from the edited populations to this graph using vg, creating a GAF (Gene Annotation Format) file representing the alignment.
A Python script is called to parse the alignment file and calculate the coverage of the homology arms compared to the reference homology arms. The output of this pipeline is a coverage table, which displays the name of each intended edit sequence along with the homology arm coverage and reference homology arm coverage.
This allows FantasticLamp to simultaneously quantify the efficacy of edits in multiple edited populations, which can be used to test novel editing methods as well as to verify current methods in experiments that involve genomic edits.
By using reads aligned to a graph instead of linear alignment, FantasticLamp can avoid reference bias, making it a useful tool for researchers in the field of genome engineering and other related fields who wish to precisely quantify genome edit states.

As novel genome editing processes are developed in the near future, we posit that the need for software tools that can capture the complexity of the edits and the relationship between them will increase.
We show a basic approach for quantifying complex genome editing results, which are difficult to reliably and simply evaluate using a single reference genome that may lead to biased estimates of edit states.
The design library CSV file used in the creation of this pipeline was unique to a specific set of editing experiments performed at Inscripta Inc. 
However, the basic format is generic to other editing experiments that users may carry out in the future, with only minimal changes to the script to adjust for difference in formatting of the design library file.
Future users will have to edit the bash script "find_coverage.sh" such that the correct columns containing the edit- and reference homology arms are extracted.
While the pipeline was tested only on singleplex sequencing data, it is possible to perform analysis using pooled sequencing data, but data availability has limited our ability to thoroughly test this. 
In summary, FantasticLamp provides a basic demonstration of the principle of using a variation graph as a reference system to avoid bias when quantifying genomic edits, and stands as a prototype for future work in this space.

# Acknowledgements
We are grateful to Inscripta Inc for providing the data from editing experiments that was used to create and test this pipeline. We thank Deanna Church for supporting our collaboration.

# References
