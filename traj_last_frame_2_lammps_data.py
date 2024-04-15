# open the input file for reading
with open('traj.lammpstrj', 'r') as f:
    # read the lines of the file into a list
    lines = f.readlines()

# get the total number of atoms from the fourth line
num_atoms = int(lines[3].split()[0])

# extract the last N+9 lines of the file
output_lines = lines[-(num_atoms+9):]

# open the output file for writing
with open('traj_last_frame.lammpstrj', 'w') as f:
    # write the extracted lines to the output file
    f.writelines(output_lines)
    f.close()

# open the input file for reading
with open('traj_last_frame.lammpstrj', 'r') as f:
    # read the lines of the file into a list
    lines = f.readlines()

# get the number of atoms
num_atoms = int(lines[3].split()[0])

# get the box bounds
xlo, xhi = map(float, lines[5].split())
ylo, yhi = map(float, lines[6].split())
zlo, zhi = map(float, lines[7].split())

# create a dictionary to map element symbols to element types
element_dict = {'C': 1, 'N': 2}

# create a list to hold the atom data
atom_data = []

# loop through the lines containing the atom data
for line in lines[9:]:
    # split the line into its fields
    fields = line.split()

    # get the atom ID, type, and element
    atom_id, atom_type, element = int(fields[0]), int(fields[1]), fields[2]

    # get the atom coordinates
    x, y, z = map(float, fields[3:6])

    # append the atom data to the list
    atom_data.append((atom_id, element_dict[element], x, y, z))

# open the output file for writing
output_file = input("Enter the output file name: ")
with open(output_file, 'w') as f:
    # write the header
    f.write('# convert by traj_last_frame_2_lammps_data.py\n\n')
    f.write('{} atoms\n'.format(num_atoms))
    f.write('1 atom types\n\n')

    # write the box bounds
    f.write('{:.8f} {:.8f} xlo xhi\n'.format(xlo, xhi))
    f.write('{:.8f} {:.8f} ylo yhi\n'.format(ylo, yhi))
    f.write('{:.8f} {:.8f} zlo zhi\n\n'.format(zlo, zhi))

    # write the Masses section
    f.write('Masses\n\n')
    f.write('1 12.0107 # C\n\n')

    # write the Atoms section
    f.write('Atoms  # atomic\n\n')
    for atom in atom_data:
        f.write('{:d} {:d} {:.6f} {:.6f} {:.6f}\n'.format(*atom))
