# Flag to indicate if we are in the desired section
# extract_lines = False

# Open the input file for reading
# with open("log.lammps",'r') as input_file:
#     # Open the output file for writing
#     with open("data.log.lammps", 'w') as output_file:
#         # Iterate over each line in the input file
#         for line in input_file:
#             # Check if the line starts with "Step Time"
#             if line.startswith("Step Time"):
#                 # Set the flag to True to start extracting lines
#                 extract_lines = True
#             # Check if the line starts with "Loop time"
#             elif line.startswith("Loop time"):
#                 # Set the flag to False to stop extracting lines
#                 extract_lines = False
#                 # Exit the loop since we've reached the end of the desired section
#                 break  

#             # If the flag is True, write the line to the output file
#             if extract_lines and not line.startswith("Step Time"):
#                 output_file.write(line)
		
		
extract_lines = False
lines_to_extract = []

with open("log.lammps", 'r') as input_file:
    for line in input_file:
        if line.startswith("Step Time"):
            # Clear previous lines only when extract_lines is True
            if extract_lines:
                lines_to_extract = []
            extract_lines = True
        elif line.startswith("Loop time"):
            extract_lines = False

        if extract_lines and not line.startswith(("Step Time", "Loop time")):
            lines_to_extract.append(line)

with open("data.log.lammps", 'w') as output_file:
    output_file.writelines(lines_to_extract)
