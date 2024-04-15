# Ask the user for inputs
def change_coordinates():

    # Ask the user for inputs
    original_file = input('Please provide the path to the original xyz file: ')
    changed_file = input('Please provide the path to the output xyz file: ')
    change_factor = float(input('Please provide the change factor: '))

    # Read the original xyz file
    with open(original_file, "r") as f:
        lines = f.readlines()

    # Get number of atoms
    num_atoms = int(lines[0])

    # Get atom data
    atom_data = [line.split() for line in lines[2:]]
    atom_data = [[d[0], float(d[1]), float(d[2]), float(d[3])] for d in atom_data]

    # Change the atom's x, y, and z coordinates according to the change_factor
    atom_data = [[d[0], change_factor*d[1], change_factor*d[2], change_factor*d[3]] for d in atom_data]

    # Change the first column of all_atoms to "C"
    for column in atom_data:
        column[0] = "C"

    # Write the new structure to file
    with open(changed_file, "w") as f:
        f.write(str(int(num_atoms)) + "\n")
        f.write("\n")
        for atom in atom_data:
            f.write("{} {:.6f} {:.6f} {:.6f}\n".format(atom[0], atom[1], atom[2], atom[3]))
    print('Operation completed successfully.')

# Run the function when the script is called from command line
if __name__ == '__main__':
    change_coordinates()
