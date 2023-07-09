The following files are strictly necessary:
DesignLibraryDetails_ODD126.withEditWindow.csv
Data_names.txt
environment.yaml
ref_and_mt.fna
simple_test.fastq.gz

Copy the following scripts into the data folder:
find_coverage.sh
compare_coverage_read_info.py

ODD126_augmented_CB39.fasta is not strictly necessary, but there will be an error message 
if the pipeline does not find it. However, the pipeline will still run correctly, as this test 
does not include reads from a vector plasmid sequence.

The pipeline needs to be run with Conda
Install deps:
conda env update --file environment.yaml

Run pipeline: (Use this exact command)
conda run -n fantastic-lamp bash find_coverage.sh

IMPORTANT: Any errors regarding "segmentation fault" and/or "core dumped" can be assumed to arise from
not exporting jemalloc correctly, normally when it is not installed or the path is incorrect.
You may need to run:
sudo apt install libjemalloc-dev if the problem persists.

Visualising graph:
To visualise the graph in BandageNG (https://github.com/asl/BandageNG), load the simple_test.og.gfa file.
To label the nodes with their coverages and other information, BandageNG requires a CSV file. To create a CSV file
from the TSV output, simply run the BASH script:
convert_to_csv.sh
in the /Test folder