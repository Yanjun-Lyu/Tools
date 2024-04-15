# Import the necessary libraries
import math

# Get user input for molar weight of particle, packing density, and volume of computational box 
molar_weight = float(input("Enter the molar weight of the particle (in g/mol): ")) 
packing_density = float(input("Enter the packing density of the particle (in g/cm^3): ")) 
natoms = float(input("Enter number of atoms to pack: ")) 

# Calculate the molecular weight of the particle
NAV = 6.02214076e23 # Avogadro's number in mol^-1
molecularMass = molar_weight / NAV 

# Convert the volume to cubic centimeters (1 Angstrom^3 = 1e-24 cm^3)
boxVol = natoms * molecularMass / packing_density  # in cm^3
boxDim = boxVol ** (1/3) * 1e8  # in angstrom

# Print the result
print(f"At density {packing_density} gcc, a cubic box size required to pack {natoms} atoms is {boxDim} angstrom")

