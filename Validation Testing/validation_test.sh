#!/bin/bash
env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'
export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2
# This script will be used to run the test in one go, tying together all the other scripts.

# Number of edits:
num_edits=101
genome_name="GCA_000146045.2_R64_genomic.fna"
# Create edits:
function create_edits() {
  # Create the edits
  python3 create_edits.py $num_edits $genome_name
}

# Create the reads
function create_reads() {
  for ((i=0; i<$num_edits; i++)); do
    wgsim -e $1 -r 0 -R 0 -X 0 -h -S 0 -d 200 -N $2 -1 200 -2 200 -S 5 genome${i}.fasta genome${i}_1.fastq genome${i}_2.fastq
  done
  # $1 is error rate (out of 1), $2 is number of (pairs of) reads
}

function zip_reads() {
  for ((i=0; i<$num_edits; i++)); do
    gzip genome${i}_1.fastq
    gzip genome${i}_2.fastq
    echo "Zipped genome${i}"
  done
}

function  create_edit_graph() {
  cat genomic_edits.csv | awk -F',' '{print ">homology_arm_"$1; print $9;}' | tr -d '\r' > homology_arms.fa
  cat genomic_edits.csv | awk -F',' '{print ">ref_homology_arm_"$1; print $8;}' | tr -d '\r' > ref_subpaths.fa
  cat ref_subpaths.fa homology_arms.fa > ref_and_hom_arms.fa
  ref=$1
  minimap2 -k 19 -w 1 -cx sr $ref ref_and_hom_arms.fa >ref_and_hom_arms.paf
  cat $ref ref_and_hom_arms.fa > ref+edits.fa
  seqwish -g ref+edits.gfa -s ref+edits.fa -p ref_and_hom_arms.paf -P
  odgi build -g ref+edits.gfa -o -  | odgi sort -i - -p sYYgs -o - | odgi chop -i - -o - -c 256 | tee ref+edits.og |odgi view -i - -g >ref+edits.og.gfa
  vg convert -x ref+edits.og.gfa >ref+edits.og.gfa.xg
  vg index -p -g ref+edits.og.gfa.gcsa -t 16 ref+edits.og.gfa.xg
  for ((i=0; i<$num_edits; i++)); do
    vg map -x ref+edits.og.gfa.xg -g ref+edits.og.gfa.gcsa -t 16 -% -f genome${i}_1.fastq.gz -f genome${i}_2.fastq.gz | pv -l >genome${i}.gaf
    python3 compare_coverage_read_info.py --gaf-path genome${i}.gaf --out-path genome${i}.tsv --og-path "ref+edits.og"
  done
}

function collect_results() {
    python3 collect_results.py $1 $num_edits
}

function summarize_coverage_results() {
  echo "Coverage Results Summary:"
  for file in *x*.csv; do
    edit_rate=$(echo $file | sed -E 's/e([0-9]+)N.*/\1/')
    total_edits=$(($(cat $file | wc -l) - 2)) # -2 for header and the negative control
    true_count=$(awk -F',' '$2=="True" {count++} END {print count}' $file)

    echo "For $file:"
    echo "  Total Edits: $total_edits"
    echo "  True Count  : $true_count"
    echo "  Success Rate: $((true_count * 100 / total_edits))%"
    echo ""
  done
}

function summarize_error_results() {
  echo "Error Results Summary:"
  for file in e*.csv; do
    edit_rate=$(echo $file | sed -E 's/e([0-9]+).*/\1/')
    total_edits=$(($(cat $file | wc -l) - 2))
    true_count=$(awk -F',' '$2=="True" {count++} END {print count}' $file)

    echo "For $file:"
    echo "  Total Edits: $total_edits"
    echo "  True Count  : $true_count"
    echo "  Success Rate: $((true_count * 100 / total_edits))%"
    echo ""
  done
}

create_edits

# Varying coverage (number of reads) with 0.00 error rate:
# 1x coverage:
create_reads 0.00 30000
zip_reads
create_edit_graph $genome_name
collect_results 1xN30K $num_edits
rm *.gz *.gaf

# 2x coverage:
create_reads 0.00 60000
zip_reads
create_edit_graph $genome_name
collect_results 2xN60K $num_edits
rm *.gz *.gaf

# 3x coverage:
create_reads 0.00 90000
zip_reads
create_edit_graph $genome_name
collect_results 3xN90K $num_edits
rm *.gz *.gaf

# 4x coverage:
create_reads 0.00 120000
zip_reads
create_edit_graph $genome_name
collect_results 4xN120K $num_edits
rm *.gz *.gaf

# 5x coverage:
create_reads 0.00 150000
zip_reads
create_edit_graph $genome_name
collect_results 5xN150K $num_edits
rm *.gz *.gaf

# 6x coverage:
create_reads 0.00 180000
zip_reads
create_edit_graph $genome_name
collect_results 6xN180K $num_edits
rm *.gz *.gaf

# 7x coverage:
create_reads 0.00 210000
zip_reads
create_edit_graph $genome_name
collect_results 7xN210K $num_edits
rm *.gz *.gaf

# 8x coverage:
create_reads 0.00 240000
zip_reads
create_edit_graph $genome_name
collect_results 8xN240K $num_edits
rm *.gz *.gaf

# 9x coverage:
create_reads 0.00 270000
zip_reads
create_edit_graph $genome_name
collect_results 9xN270K $num_edits
rm *.gz *.gaf
##
### 10x coverage:
create_reads 0.00 300000
zip_reads
create_edit_graph $genome_name
collect_results 10xN300K $num_edits
rm *.gz *.gaf

# Varying base error rate (with 10x coverage):
# 0.000
create_reads 0.00 300000
zip_reads
create_edit_graph $genome_name
collect_results e000 $num_edits
rm *.gz *.gaf
#
#0.001
create_reads 0.001 300000
zip_reads
create_edit_graph $genome_name
collect_results e001 $num_edits
rm *.gz *.gaf

#0.005
create_reads 0.005 300000
zip_reads
create_edit_graph $genome_name
collect_results e005 $num_edits
rm *.gz *.gaf

#0.01
create_reads 0.01 300000
zip_reads
create_edit_graph $genome_name
collect_results e01 $num_edits
rm *.gz *.gaf

#0.02
create_reads 0.02 300000
zip_reads
create_edit_graph $genome_name
collect_results e02 $num_edits
rm *.gz *.gaf

#0.05
create_reads 0.05 300000
zip_reads
create_edit_graph $genome_name
collect_results e05 $num_edits
rm *.gz *.gaf

#0.1
create_reads 0.1 300000
zip_reads
create_edit_graph $genome_name
collect_results e1 $num_edits
rm *.gz *.gaf

summarize_coverage_results
summarize_error_results

echo "Done!"
