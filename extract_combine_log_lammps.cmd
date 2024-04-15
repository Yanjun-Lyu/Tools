#!/bin/bash

# Set the working directory as the directory where the script is executed
working_dir="$PWD"

# Check if the maximum number of runs is provided as a command-line argument
if [ $# -eq 0 ]; then
    echo "Please provide the maximum number of runs as a command-line argument"
    exit 1
fi

# Get the maximum number of runs from the command-line argument
max_run=$1

# Initialize the combined file
echo -n > data_1_${max_run}.log.lammps

# Loop over all the directories from run_1 to run_max
for ((i=1; i<=$max_run; i++))
do
        # Change to the directory
        cd "$working_dir/run_${i}"

        # Extract data from log.lammps into data.log.lammps
        extract_data_log_lammps='python3 /usr/workspace/wsa/lyu1/utils/extract_data_from_log_lammps.py'
        ${extract_data_log_lammps}

        # Append the data to the combined file
        cat data.log.lammps >> "$working_dir/data_1_${max_run}.log.lammps"

        # Move back to the parent directory
        cd "$working_dir"
done

