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
  //////////////////////////////////////////////////////
  // Checking input errors
  //////////////////////////////////////////////////////

  // Exit execution if input argument is incorrect
  if (argc != 4)
    {
      cerr << "Usage: ./<this_file> <traj_file> <binWidth> <binning direction (1 for x, 2 for z)> \n";   // cerr will not be flushed.
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
  vector<double> centerBinCoords; // x (or z) coordinate of each bin index in one frame

  vector<vector<double>> denBinAllFrame; // Record all coordination number profile for calculation of standard deviation
  denBinAllFrame.clear();

  vector<double> avgCenterBinCoords;    // Average x (or z) coordinate of each bin index in all frames
  vector<double> avgDenBin; // total coordination number of atom in the bin in all frames

  cout << endl;
  cout << "# Input parameters:" << endl;
  cout <<                             "# Trajectory file:   " << TRAJ_FILE << endl;  
  cout << fixed << setprecision(4) << "# Bin width:         " << binWidth << endl; 
  cout <<                             "# Binning direction: " << binDir << endl; 


  ////////////////////////////////////////////////////////////
  // Setting output files
  ////////////////////////////////////////////////////////////

  // For writing density profile into a file
  const string DEN_PROFILE = "den_profile_" + TRAJ_FILE + ".dat";
  string rm  = "rm -f " + DEN_PROFILE;
  system (rm.c_str());
  ofstream denProfile(DEN_PROFILE, ios::app);

  // For writing density profile average (all frames) into a file
  const string DEN_PROFILE_AVG = "den_profile_avg_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + DEN_PROFILE_AVG;
  system(rm.c_str());
  ofstream denProfileAvg(DEN_PROFILE_AVG, ios::app);


  ////////////////////////////////////////////////////////////
  // Reading trajectory file
  ////////////////////////////////////////////////////////////

  int nframes = 0;
  bool calledBefore = false;

  cout << endl;

  double actualBinWidth;  // for the last bin which is smaller than the regular one


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
        centerBinCoords.resize(nbins, 0.0);
        denBin.resize(nbins, 0.0);

        avgDenBin.resize(nbins, 0.0);
        avgCenterBinCoords.resize(nbins, 0.0);
 
        calledBefore = true;
      }

      // Wrap coordinates
      traj.pbc_wrap(); 

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
      // output density profile into dat file
      denProfile << "# timestep: " << traj.tstep << endl; 
      denProfile << "# center_x_coord  avg_mass_den (gcc)" << endl; 
      for (size_t i = 0; i < nbins; i++)
      {

        if (i == nbins - 1)
        {
          centerBinCoords[i] = (traj.boxdims.x - i * binWidth) * 0.5 + i * binWidth + traj.bounds[0].x;
        }
        else
        {
          centerBinCoords[i] = traj.bounds[0].x + (i + 0.5) * binWidth;
        }

        actualBinWidth = min(binWidth, traj.boxdims.x - (i - 1) * binWidth);

        denBin[i] = natomsBin[i] * MASS_C_ATOM / (traj.boxdims.y * traj.boxdims.z * actualBinWidth * pow(A2CM, 3));       // Unit: gcc
                  
        avgCenterBinCoords[i] += centerBinCoords[i];
        
      }

      for (size_t i = 1; i < nbins - 1; i++)
      {
        denProfile << centerBinCoords[i] << "  " << denBin[i] << endl; 
      }

      // Record density profiles of the frame for calculation of average of bin among all frames
      denBinAllFrame.push_back(denBin);

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

        natomsBin.resize(nbins,0);
        centerBinCoords.resize(nbins, 0.0);
        denBin.resize(nbins, 0.0);

        avgDenBin.resize(nbins, 0.0);
        avgCenterBinCoords.resize(nbins, 0.0);
 
        calledBefore = true;
      }

      // Wrap coordinates
      traj.pbc_wrap(); 

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
      // output density profile into dat file
      denProfile << "# timestep: " << traj.tstep << endl; 
      denProfile << "# center_z_coord  avg_mass_den (gcc)" << endl; 
      for (size_t i = 0; i < nbins; i++)
      {

        if (i == nbins - 1)
        {
          centerBinCoords[i] = (traj.boxdims.z - i * binWidth) * 0.5 + i * binWidth + traj.bounds[0].z;
        }
        else
        {
          centerBinCoords[i] = traj.bounds[0].z + (i + 0.5) * binWidth;
        }

        actualBinWidth = min(binWidth, traj.boxdims.z - (i - 1) * binWidth);

        denBin[i] = natomsBin[i] * MASS_C_ATOM / (traj.boxdims.x * traj.boxdims.y * actualBinWidth * pow(A2CM, 3));       // Unit: gcc
                  
        avgCenterBinCoords[i] += centerBinCoords[i];
        
      }
      
      for (size_t i = 1; i < nbins - 1; i++)
      {
        denProfile << centerBinCoords[i] << "  " << denBin[i] << endl; 
      }

      // Record density profiles of the frame for calculation of average of bin among all frames
      denBinAllFrame.push_back(denBin);

      nframes++;

      cout << endl;
    }
    break;

  }


  // Calculate the average center coordinate of each bin among all frames
  // Calculate average of the quantities of each bin among all frames
  // Note: dimension of the "all-frame" vectors: (nrow = nframes, ncol = nbins) 
  vector<double> BinDenBinOfAllFrame; // column vector of the all-frame vector
  BinDenBinOfAllFrame.resize(nframes, 0.0); 

  for (size_t j = 0; j < nbins; j++)
  {
    avgCenterBinCoords[j] /= nframes; 
    for (size_t i = 0; i < nframes; i++)
    {
      BinDenBinOfAllFrame[i] = denBinAllFrame[i][j];
    }

    avgDenBin[j] = mean(BinDenBinOfAllFrame);
  }
  
  // Output average density profile of all frame
  denProfileAvg << "# number of timesteps averaged: " << nframes << endl; 
  denProfileAvg << "# avg_center_coord  avg_coord_num_all_frames" << endl;

  for (size_t i = 1; i < nbins - 2; i++) // skip the last 2 bins (adjust as needed)
  {
    denProfileAvg << fixed << setprecision(4) << avgCenterBinCoords[i] << "  " << avgDenBin[i] << endl;
  }

  cout << endl;
  cout << "# " << nframes  << " frames are read and processed." << endl;
  cout << endl;

  denProfile.close();
  cout << "# Mass density profile written to " << DEN_PROFILE << endl;

  denProfileAvg.close();
  cout << "# Average mass density profile written to " << DEN_PROFILE_AVG << endl;
  
  return 0;
}

