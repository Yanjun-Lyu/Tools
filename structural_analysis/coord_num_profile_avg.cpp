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

  unsigned int binDir = stoi(argv[4]); // Binning direction (x, y or z)
  
  if (binDir != 1 && binDir != 2)
  {
    cerr << "Invalid binning direction. It should either be 1 (x) or 2 (z)." << endl;
    return 1;
  }  

  vector<int> natomsBin;        // Number of atoms in the bin of interest
  vector<double> coordNumBin;   // total coordination number of atom in the bin in one frame 
  vector<double> centerXCoords;  // Average x coordinate of each bin index in one frame

  vector<double> avgCenterXCoords;    // Average x coordinate of each bin index in all frames
  vector<double> coordNumBinAllFrame;   // total coordination number of atom in the bin in all frames

  ////////////////////////////////////////////////////////////
  // Setting output files
  ////////////////////////////////////////////////////////////

  // For writing coordination number profile into a file
  const string COORD_NUM_PROFILE = "coord_num_profile_" + TRAJ_FILE + ".dat";
  string rm = "rm -f " + COORD_NUM_PROFILE;
  system(rm.c_str());
  ofstream coordNumProfile(COORD_NUM_PROFILE, ios::app);

  // For writing coordination number profile average into a file
  const string COORD_NUM_PROFILE_AVG = "coord_num_profile_avg_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + COORD_NUM_PROFILE_AVG;
  system(rm.c_str());
  ofstream coordNumProfileAvg(COORD_NUM_PROFILE_AVG, ios::app);


  ////////////////////////////////////////////////////////////
  // Reading and processing trajectory file
  ////////////////////////////////////////////////////////////

  int nframes = 0;
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

          natomsBin.resize(nbins, 0);
          centerXCoords.resize(nbins, 0.0);
          coordNumBin.resize(nbins, 0.0);
          avgCenterXCoords.resize(nbins, 0.0);
          coordNumBinAllFrame.resize(nbins, 0.0);
          
          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap(); 

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
        // Calculate center X coordinate number of each bin in the frame
        // Populate coordNumBinAllFrame and avgCenterXCoords for all frames
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumBin[i] /= natomsBin[i];
          centerXCoords[i] = traj.bounds[0].x + (i + 0.5) * binWidth; 
          
          coordNumBinAllFrame[i] += coordNumBin[i];
          avgCenterXCoords[i]  += centerXCoords[i]; 
        }
        
        // output coordination number profile into dat file
        coordNumProfile << "# timestep: " << traj.tstep << endl; 
        coordNumProfile << "# center_x_coord  avg_coord_num" << endl; 
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumProfile << fixed << setprecision(8) << centerXCoords[i] << "  " << coordNumBin[i] << endl;
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
          nbins = traj.boxdims.z / binWidth;
          cout << "# number of bins: " << nbins << endl;

          natomsBin.resize(nbins, 0);
          centerXCoords.resize(nbins, 0.0);
          coordNumBin.resize(nbins, 0.0);
          avgCenterXCoords.resize(nbins, 0.0);
          coordNumBinAllFrame.resize(nbins, 0.0);
          
          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap(); 

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
        // Calculate center X coordinate number of each bin in the frame
        // Populate coordNumBinAllFrame and avgCenterXCoords for all frames
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumBin[i] /= natomsBin[i];
          centerXCoords[i] = traj.bounds[0].z + (i + 0.5) * binWidth; 
          
          coordNumBinAllFrame[i] += coordNumBin[i];
          avgCenterXCoords[i]  += centerXCoords[i]; 
        }
        
        // output coordination number profile into dat file
        coordNumProfile << "# timestep: " << traj.tstep << endl; 
        coordNumProfile << "# center_z_coord  avg_coord_num" << endl; 
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumProfile << fixed << setprecision(8) << centerXCoords[i] << "  " << coordNumBin[i] << endl;
        }

        nframes++;

        cout << endl;
      }
      break;
  }

  // Calculate average coordination number and center x coordinates among all frames
  for (size_t i = 0; i < nbins; i++)
  {
    coordNumBinAllFrame[i] /= nframes;
    avgCenterXCoords[i]    /= nframes; 
  }

  // Output average coordination number profile into dat file
  coordNumProfileAvg << "# number of timesteps averaged: " << nframes << endl; 
  coordNumProfileAvg << "# avg_center_coord  avg_coord_num_all_frames" << endl; 
  for (size_t i = 0; i < nbins; i++)
  {
    coordNumProfileAvg << fixed << setprecision(8) << avgCenterXCoords[i] << "  " << coordNumBinAllFrame[i] << endl;
  }

  cout << endl;
  cout << "# " << nframes  << " frames are read" << endl;
  cout << endl;

  coordNumProfile.close();
  cout << "# Coordination number profile written to " << COORD_NUM_PROFILE << endl;

  coordNumProfileAvg.close();
  cout << "# Average coordination number profile written to " << COORD_NUM_PROFILE_AVG << endl;

  return 0;
}
