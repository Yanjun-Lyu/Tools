#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "helpers.hpp"
#include "trajectory.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////
  // Checking input errors
  ////////////////////////////////////////////////////////////

  // Exit execution if input argument is incorrect
  if (argc != 5)
    {
      cerr << "Usage: ./<this_file> <traj_file> <cutoff> <binWidth> <1 (bin in x) or 2 (bin in z)> \n";   // cerr will not be flushed.
	    return 1;
    }
  
  unsigned int binDir = stoi(argv[4]); // Binning direction (x, y or z)
  
  if (binDir != 1 && binDir != 2)
  {
    cerr << "Invalid binning direction. It should either be 1 (x) or 2 (z)." << endl;
    return 1;
  }
  
  const string TRAJ_FILE = argv[1];
  Trajectory traj(TRAJ_FILE);

  ////////////////////////////////////////////////////////////
  // Setting parameters
  ////////////////////////////////////////////////////////////

  // Bins for coordination number profile
  int nbins;                           // Number of bins
  int binIndex;                        // Bin index of an atom of interest
  double cutoff = stod(argv[2]);       // Cutoff distance for coord num
  double binWidth = stod(argv[3]);     // Bin width for coord num profile (Ang)
  
  vector<int> natomsBin;        // Number of atoms in the bin of interest
  vector<double> coordNumBin;   // total coordination number of atom in the bin


  ////////////////////////////////////////////////////////////
  // Setting output files
  ////////////////////////////////////////////////////////////

  // For writing wrapped coordinates into a new lammpstrj file
  const string WRAPPED_TRAJ_FILE = "wrapped_" + TRAJ_FILE;
  string rm = "rm -f " + WRAPPED_TRAJ_FILE; 
  system(rm.c_str());
  ofstream wrappedCoords(WRAPPED_TRAJ_FILE, ios::app);

  // For writing coordination number profile into a file
  const string COORD_NUM_PROFILE = "coord_num_profile_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + COORD_NUM_PROFILE; 
  system(rm.c_str());
  ofstream coordNumProfile(COORD_NUM_PROFILE, ios::app);


  ////////////////////////////////////////////////////////////
  // Reading trajectory file
  ////////////////////////////////////////////////////////////

  int nframes = 1;
  bool calledBefore = false;

  cout << endl;

  switch(binDir)
  {
    // For binning in x direction
    case 1:
      cout << "# Binning in x direction" << endl;
      cout << endl;

      while(traj.read_frame()) // Read frames until the end of file
      {
        // Calculate number of bins and initialize vectors
        if (!calledBefore)
        {
          nbins = traj.boxdims.x / binWidth;
          cout << "# number of bins: " << nbins << endl;

          natomsBin.resize(nbins,0);
          coordNumBin.resize(nbins, 0.0);
          
          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap(); 

        // output wrapped coordinates in a lammpstrj file
        wrappedCoords << "ITEM: TIMESTEP" << endl;
        wrappedCoords << traj.tstep << endl; 
        wrappedCoords << "ITEM: NUMBER OF ATOMS" << endl;
        wrappedCoords << traj.natoms << endl;
        wrappedCoords << "ITEM: BOX BOUNDS pp pp pp" << endl;
        wrappedCoords << traj.bounds[0].x << " " << traj.bounds[1].x << endl;
        wrappedCoords << traj.bounds[0].y << " " << traj.bounds[1].y << endl;
        wrappedCoords << traj.bounds[0].z << " " << traj.bounds[1].z << endl;
        wrappedCoords << "ITEM: ATOMS id type element xu yu zu" << endl; 
        for (size_t i = 0; i < traj.natoms; i++)
        {
          wrappedCoords << fixed << setprecision(8) << i + 1 << " 1 C " << traj.coords[i].x << " " << traj.coords[i].y << " " << traj.coords[i].z << endl;
        }

        // Get coordination number of each atom (rate-determining step)
        traj.get_coord_num(cutoff); 
      
        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        coordNumBin.assign(coordNumBin.size(), 0);

        // Get and populate coordination number in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].x - traj.bounds[0].x) / binWidth;
          natomsBin[binIndex]++;
          coordNumBin[binIndex] += traj.coordNums[i];
        }

        // Calculate average coordination number in each bin
        // Calculate average mass density of each bin
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumBin[i] /= natomsBin[i];
        }
        
        // output coordination number profile into dat file
        coordNumProfile << "# timestep: " << traj.tstep << endl; 
        coordNumProfile << "# center_x_coord  avg_coord_num" << endl; 
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumProfile << fixed << setprecision(8) << traj.bounds[0].x + (i + 0.5) * binWidth << "  " << coordNumBin[i] << endl;
        }

        nframes++;

        cout << endl;
      }
      break;

    // For in binning in z direction
    case 2:
      cout << "Binning in z direction" << endl;
      cout << endl;

      while(traj.read_frame()) // Read frames until the end of file
      {
        // Wrap coordinates
        traj.pbc_wrap();

        // output wrapped coordinates in a lammpstrj file
        wrappedCoords << "ITEM: TIMESTEP" << endl;
        wrappedCoords << traj.tstep << endl; 
        wrappedCoords << "ITEM: NUMBER OF ATOMS" << endl;
        wrappedCoords << traj.natoms << endl;
        wrappedCoords << "ITEM: BOX BOUNDS pp pp pp" << endl;
        wrappedCoords << traj.bounds[0].x << " " << traj.bounds[1].x << endl;
        wrappedCoords << traj.bounds[0].y << " " << traj.bounds[1].y << endl;
        wrappedCoords << traj.bounds[0].z << " " << traj.bounds[1].z << endl;
        wrappedCoords << "ITEM: ATOMS id type element xu yu zu" << endl; 
        for (size_t i = 0; i < traj.natoms; i++)
        {
          wrappedCoords << fixed << setprecision(8) << i + 1 << " 1 C " << traj.coords[i].x << " " << traj.coords[i].y << " " << traj.coords[i].z << endl;
        }

        // Calculate number of bins and initialize vectors
        if (!calledBefore)
        {
          nbins = traj.boxdims.z / binWidth;
          cout << "# number of bins: " << nbins << endl;

          natomsBin.resize(nbins,0);
          coordNumBin.resize(nbins, 0.0);
          
          calledBefore = true;
        }

        // Get coordination number of each atom (rate-determining step)
        traj.get_coord_num(cutoff); 
      
        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        coordNumBin.assign(coordNumBin.size(), 0);

        // Get and populate coordination number in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].z - traj.bounds[0].z) / binWidth;
          natomsBin[binIndex]++;
          coordNumBin[binIndex] += traj.coordNums[i];
        }

        // Calculate average coordination number in each bin
        // Calculate average mass density of each bin
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumBin[i] /= natomsBin[i];
        }
        
        // output coordination number profile into dat file
        coordNumProfile << "# timestep: " << traj.tstep << endl; 
        coordNumProfile << "# center_z_coord  avg_coord_num" << endl; 
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumProfile << fixed << setprecision(8) << traj.bounds[0].z + (i + 0.5) * binWidth << "  " << coordNumBin[i] << endl;
        }

        nframes++;

        cout << endl;
      }
      break;
  }
  
  cout << endl;
  cout << "# " << nframes - 1  << " frames are read" << endl;
  cout << endl;

  // Close file and confirm writing.
  wrappedCoords.close(); 
  cout << "# Wrapped trajectory written to: " << WRAPPED_TRAJ_FILE << endl;
  cout << endl;

  coordNumProfile.close();
  cout << "# Coordination number profile written to " << COORD_NUM_PROFILE << endl;
  
  return 0;
}









  // // test get_mass_den
  // vector<xyz> usr_sp_bounds;
  // usr_sp_bounds.push_back({0.0, 0.0, 0.0});           // (xlo, ylo, zlo)
  // usr_sp_bounds.push_back({11.471, 11.471, 11.471});     // (xhi, yhi, zhi)
  // double den = traj.get_mass_den(usr_sp_bounds);
  // cout << "Mass density: " << den << " gcc" << endl;
 
  // // test read_frame at the last line of file
  // cout << "Atom 250: " << traj.coords[249].x << " " << traj.coords[249].y << " " << traj.coords[249].z << endl;

