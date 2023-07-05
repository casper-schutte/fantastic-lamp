#!/bin/bash

input_file="simple_test.tsv"
output_file="simple_test.csv"

sed 's/\t/,/g' "$input_file" > "$output_file"


