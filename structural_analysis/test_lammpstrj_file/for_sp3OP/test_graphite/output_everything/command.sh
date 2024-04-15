#!/bin/bash

# Loop through numbers 0 to 12
for (( i=0; i<=12; i++ ))
do
    # Define output file name based on the pattern
    output_file="atomsInBin_$i.out"

    # Use awk to match lines ending with "is in Bin i"
    awk -v pattern="is in Bin $i\$" '$0 ~ pattern { print > "'"$output_file"'" }' OP.out
done



# parsing bin 0

awk '{ print $2 ":" }' atomsInBin_0.out > atomsInBin_0.out.pattern

awk 'FNR==NR { patterns[$0]; next } { for (pattern in patterns) { if ($0 ~ pattern) print } }' atomsInBin_0.out.pattern sp3OP.out > atomsInBin_0.out.sp3OP
