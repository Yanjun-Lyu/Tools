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
echo -n > traj_1_${max_run}.lammpstrj

# Loop over all the directories from run_1 to run_max
for ((i=1; i<=$max_run; i++))
do
        # Change to the directory
        cd "$working_dir/run_${i}"

        # Append the data to the combined file
        cat traj.lammpstrj >> "$working_dir/traj_1_${max_run}.lammpstrj"

        # Move back to the parent directory
        cd "$working_dir"
done

