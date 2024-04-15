#!/bin/bash

# Usage: ./avg_stddev.sh <data_file> <column_number> <num_rows_from_end> <output_file>

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <data_file> <column_number> <num_rows_from_end> <output_file>"
  read -p "Enter data file, column number, num rows from end, and output file (separated by spaces): " data_file column_number num_rows_from_end output_file
fi

data_file=$1
column_number=$2
num_rows_from_end=$3
output_file=$4

# Calculate the average using awk
average=$(tail -n "$num_rows_from_end" "$data_file" | awk -v col="$column_number" '{ sum += $col; count++ } END { if (count > 1) printf "    %6.4f", sum/count; else print "NAN" }')

# Calculate the standard deviation using awk
standard_deviation=$(tail -n "$num_rows_from_end" "$data_file" | awk -v col="$column_number" -v avg="$average" '{ sum += ($col - avg) ** 2; count++ } END { if (count > 1) printf "    %6.4f", sqrt(sum/(count-1)); else print "NAN" }')

# Write results to the output file
echo "# Average  Standard Deviation" > "$output_file"
echo "$average $standard_deviation" >> "$output_file"

echo "Results saved to $output_file"

