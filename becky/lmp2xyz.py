import sys
import subprocess
from math import floor

# Usage is: python <this_script> <units> <lammpstrj file> <#atoms> <optional: frame skip>
#
# units are either: "REAL" or "METAL" or "NATIVE" ; NATIVE preserves input file units
#
# Expects the following custom format: ATOMS id type element xu yu zu <anything else>
# Prints only coordinates, and forces IF line has fx fy fz entries... 
# Converts from REAL or METAL units to chimes_lsq and chimes_md units 
#
#
# WARNING: Assumes(and expects) an orthorhombic box
# WARNING: Ancillary support for units METAL
#
# Takes about 1 minute per 1000 frames for 12500 atoms

UNITS    = sys.argv[1]
INFILE   = sys.argv[2]
ATOMS    = int(sys.argv[3])
OUTFILE  = INFILE + ".xyz"
IFSTREAM = open(INFILE,  'r')
OFSTREAM = open(OUTFILE, 'w')

SKIP = 1

if len(sys.argv) == 5:
        SKIP = int(sys.argv[4])

print "Processing file:          " + INFILE
print "Writing output to:        " + OUTFILE
print "Expected atoms:           " + `ATOMS`

# Count the number of lines in the lammps file, convert to a number of frames

FRAMES        = 0
GET_LINES     = subprocess.Popen(["wc", "-l", INFILE], stdout=subprocess.PIPE)
(FRAMES, err) = GET_LINES.communicate()

LINES = FRAMES.split()[0]

print "Counted lines:            " + LINES

FRAMES        = floor(int(FRAMES.split()[0])/(ATOMS+9)) # -1

if FRAMES*(ATOMS+9) > LINES:
        FRAMES -= 1

print "Counted frames:           " + `FRAMES`
print "Printing every nth frame: " + `SKIP`

# Process the file

BOX_X = 0.0
BOX_Y = 0.0
BOX_Z = 0.0

CONVERSION = 1.0/1.88973/627.509551  # For lammps units real

if UNITS == "METAL":
        CONVERSION = 1.0/1.88973/27.211399

if UNITS == "NATIVE":
        CONVERSION = 1.0

PRINTED = 0

for i in xrange(int(FRAMES)):


        if (i+1)%SKIP == 0:
                OFSTREAM.write(`ATOMS` + '\n')
                PRINTED += 1

        # Ignore un-neded headers

        IFSTREAM.readline() # ITEM: TIMESTEP
        IFSTREAM.readline() # <timestep>
        IFSTREAM.readline() # ITEM: NUMBER OF ATOMS
        IFSTREAM.readline() # <atoms>
        IFSTREAM.readline() # ITEM: BOX BOUNDS xy xz yz pp pp pp

        # Save/print the box lengths

        TMP_X = IFSTREAM.readline().split()
        TMP_Y = IFSTREAM.readline().split()
        TMP_Z = IFSTREAM.readline().split()

        BOX_X = str( float(TMP_X[1]) - float(TMP_X[0]) )
        BOX_Y = str( float(TMP_Y[1]) - float(TMP_Y[0]) )
        BOX_Z = str( float(TMP_Z[1]) - float(TMP_Z[0]) )

        if (i+1)%SKIP == 0:
                OFSTREAM.write(BOX_X + " " + BOX_Y + " " + BOX_Z + '\n')

        LINE = IFSTREAM.readline() # ITEM: ATOMS id type element xu yu zu
        LINE = LINE.split()

        keys = ["element", "xu", "yu", "zu", "fx", "fy", "fz"]
        locs = [None]*len(keys)

        for j in xrange(len(keys)):

                try:
                        locs[j] = LINE.index(keys[j])-2
                except:
                        locs[j] = -1
                        #print "Warning: Property keys",keys[j],"not found ... ignoring"


        PRINT_F = True

        if (locs[0] == -1):
                print "Error: id missing from file"
                exit()

        if (locs[1] == -1) or (locs[2] == -1) or (locs[3] == -1):
                print "Error: xu, yu, or zu missing from file"
                exit()

        if (locs[4] == -1) or (locs[5] == -1) or (locs[6] == -1):

                if PRINT_F and (i == 0):
                        print "Warning: fx, fy, or fz missing from file"
                PRINT_F = False

                if i == 0:
                        print "Will not print forces"

        # Write the coordinates

        for j in xrange(ATOMS):

                LINE = IFSTREAM.readline().split()

                if (i+1)%SKIP == 0:

                        TMP = LINE[locs[0]] + " " + LINE[locs[1]] + " " + LINE[locs[2]] + " " + LINE[locs[3]]

                        if PRINT_F:
                                TMP += " " + `float(LINE[locs[4]])*CONVERSION` + " " + `float(LINE[locs[5]])*CONVERSION` + " " + `float(LINE[locs[6]])*CONVERSION`

                        OFSTREAM.write(TMP + '\n')

OFSTREAM.close()
IFSTREAM .close()

print "Printed frames:           " + `PRINTED`
