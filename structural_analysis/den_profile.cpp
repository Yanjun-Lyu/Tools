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
  if (argc != 3)
    {
      cerr << "Usage: ./<this_file> <traj_file> <binWidth> <1 (bin in x) or 2 (bin in z) \n";   // cerr will not be flushed.
	    return 1;
    }

  unsigned int binDir = stoi(argv[3]); // Binning direction (x, y or z)
  
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

  // Constants
  const double MASS_C_ATOM = 1.9944733e-23;         // A carbon atom mass in g
  const double A2CM = 1e-8;                         // Angstrom to centimeter

  // Bins for coordination number profile
  int nbins;                           // Number of bins
  int binIndex;                        // Bin index of an atom of interest
  double binWidth = stod(argv[2]);     // Bin width for coord num profile (Ang)
  
  vector<int> natomsBin;        // Number of atoms in the bin of interest
  vector<double> denBin;        // mass density of a bin


  ////////////////////////////////////////////////////////////
  // Setting output files
  ////////////////////////////////////////////////////////////

  // For writing wrapped coordinates into a new lammpstrj file
  const string WRAPPED_TRAJ_FILE = "wrapped_" + TRAJ_FILE;
  string rm = "rm -f " + WRAPPED_TRAJ_FILE; 
  system(rm.c_str());
  ofstream wrappedCoords(WRAPPED_TRAJ_FILE, ios::app);

  // For writing density profile into a file
  const string DEN_PROFILE = "density_profile_" + TRAJ_FILE + ".dat";
  rm  = "rm -f " + DEN_PROFILE;
  system (rm.c_str());
  ofstream denProfile(DEN_PROFILE, ios::app);


  ////////////////////////////////////////////////////////////
  // Reading trajectory file
  ////////////////////////////////////////////////////////////

  int nframes = 1;
  bool calledBefore = false;

  cout << endl;


  // For binning in x direction
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
        denBin.resize(nbins, 0.0);
        
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
        wrappedCoords << i + 1 << " 1 C " << traj.coords[i].x << " " << traj.coords[i].y << " " << traj.coords[i].z << endl;
      }

      // Initialize vars
      natomsBin.assign(natomsBin.size(), 0);
      denBin.assign(denBin.size(), 0.0);

      // Populate atom number in each bin
      for (size_t i = 0; i < traj.natoms; i++)
      {
        binIndex = (traj.coords[i].x - traj.bounds[0].x) / binWidth;
        natomsBin[binIndex]++;
      }

      // Calculate average mass density of each bin
      for (size_t i = 0; i < nbins; i++)
      {
        denBin[i] = natomsBin[i] * MASS_C_ATOM / (traj.boxdims.y * traj.boxdims.z * binWidth * pow(A2CM, 3));       // Unit: gcc
      }
      
      // output density profile into dat file
      denProfile << "# timestep: " << traj.tstep << endl; 
      denProfile << "# center_x_coord  avg_mass_den (gcc)" << endl; 
      for (size_t i = 0; i < nbins; i++)
      {
        denProfile << traj.bounds[0].x + (i - 0.5) * binWidth << "  " << denBin[i] << endl;
      }

      nframes++;

      cout << endl;
    }
    break;


    // For binning in z direction
    case 2:
    cout << "# Binning in z direction" << endl;
    cout << endl;

    while(traj.read_frame()) // Read frames until the end of file
    {
      // Calculate number of bins and initialize vectors
      if (!calledBefore)
      {
        nbins = traj.boxdims.x / binWidth;
        cout << "# number of bins: " << nbins << endl;

        natomsBin.resize(nbins,0);
        denBin.resize(nbins, 0.0);
        
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
        wrappedCoords << i + 1 << " 1 C " << traj.coords[i].x << " " << traj.coords[i].y << " " << traj.coords[i].z << endl;
      }

      // Initialize vars
      natomsBin.assign(natomsBin.size(), 0);
      denBin.assign(denBin.size(), 0.0);

      // Populate atom number in each bin
      for (size_t i = 0; i < traj.natoms; i++)
      {
        binIndex = (traj.coords[i].z - traj.bounds[0].z) / binWidth;
        natomsBin[binIndex]++;
      }

      // Calculate average mass density of each bin
      for (size_t i = 0; i < nbins; i++)
      {
        denBin[i] = natomsBin[i] * MASS_C_ATOM / (traj.boxdims.x * traj.boxdims.y * binWidth * pow(A2CM, 3));       // Unit: gcc
      }
      
      // output density profile into dat file
      denProfile << "# timestep: " << traj.tstep << endl; 
      denProfile << "# center_z_coord  avg_mass_den (gcc)" << endl; 
      for (size_t i = 0; i < nbins; i++)
      {
        denProfile << traj.bounds[0].z + (i - 0.5) * binWidth << "  " << denBin[i] << endl;
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

  denProfile.close();
  cout << "# Mass density profile written to " << DEN_PROFILE << endl;
  
  return 0;
}

