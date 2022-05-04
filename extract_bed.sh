#!/bin/bash

# Got the H-arm names from my the csv file where they are just numbered. Added the "homology_arm_" 
# bit. 


odgi untangle -R hom_names.txt -i yeast+edits.og -t 4 > test.bed

cat test.bed | awk -F'/t' '{print $1,$2,$3,"ref_"$4}' > test2.bed 

cat test2.bed | awk '{if (($1 ~ /I/) || ($1 ~ /V/) || ($1 ~ /X/)) print $1,$2,$3,$4}' OFS="\t" > refs.bed

# Next time, use awk '!($1 ~ /^homology/)' 

# Tried using odgi inject but there are duplicated names in the refs.bed file
