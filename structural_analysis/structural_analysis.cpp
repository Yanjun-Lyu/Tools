/*
This is the integrated structural analysis code for carbon including:
1. coordination number profile (with average and standard deviation)
2. sp3 and sp2 order parameters (with average and standard deviation)
*/

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

int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////
  // Checking input errors
  ////////////////////////////////////////////////////////////

  // Exit execution if input argument is incorrect
  if (argc != 7)
    {
      cerr << "Usage: ./<this_file> <traj_file> <neighlist_cutoff> <profile_bin_width> <sp3OP_delta_theta> <sp2OP_delta_theta> <binning_direction (1 for x, 2 for z)>\n";   // cerr will not be flushed.
	    return 1;
    }


  ////////////////////////////////////////////////////////////
  // Setting parameters
  ////////////////////////////////////////////////////////////

  // Trajectory file
  const string TRAJ_FILE = argv[1];  // Only 4 atoms in the test traj file
  Trajectory traj(TRAJ_FILE); // Create the traj object

  // Cutoff distance for building neighbor list and calculating sp2OP
  double cutoff = stod(argv[2]);

  // Bins for orderparameter profile
  int nbins;                           // Number of bins
  int binIndex;                        // Bin index of an atom of interest
  double binWidth = stod(argv[3]);     // Bin width for sp2OP profile (Ang)

  // Tolerance angle theta for order parameters
  double sp3OP_del_theta = stod(argv[4]);
  double sp2OP_del_theta = stod(argv[5]);

  // Binning direction (x, y or z)
  unsigned int binDir = stoi(argv[6]); 

  vector<int> natomsBin;        // Number of atoms in the bin of interest
  vector<double> coordNumBin;   // total coordination number of atom in the bin in one frame 
  vector<double> sp2OPBin;      // total sp2OP of atom in the bin (populated)
  vector<double> sp3OPBin;      // total sp3OP of atom in the bin (populated)
  vector<double> centerBinCoords; // x (or z) coordinate of each bin index in one frame
  
  vector<vector<double>> coordNumBinAllFrame; // Record all coordination number profile for calculation of standard deviation
  coordNumBinAllFrame.clear();
  vector<vector<double>> sp2OPBinAllFrame; // Record all sp2OP profile for calculation of standard deviation 
  sp2OPBinAllFrame.clear();
  vector<vector<double>> sp3OPBinAllFrame; // Record all sp3OP profile for calculation of standard deviation
  sp3OPBinAllFrame.clear();

  vector<double> avgCenterBinCoords;    // Average x (or z) coordinate of each bin index in all frames
  vector<double> avgCoordNumBin; // total coordination number of atom in the bin in all frames
  vector<double> avgSp2OPBin;    // total sp2OP of atoms in the bin in all frames
  vector<double> avgSp3OPBin;    // total sp3OP of atoms in the bin in all frames

  vector<double> stdCoordNumBin; // standard deviation of coordination number in a bin in all frames
  vector<double> stdSp2OPBin; // standard deviation of sp2OP in a bin in all frames
  vector<double> stdSp3OPBin; // standard deviation of sp3OP in a bin in all frames

  cout << endl;
  cout << "# Input parameters:" << endl;
  cout <<                             "# Trajectory file:   " << TRAJ_FILE << endl;  
  cout << fixed << setprecision(4) << "# neighlist cutoff:  " << cutoff << endl; 
  cout << fixed << setprecision(4) << "# sp3OP del_theta:   " << sp3OP_del_theta << endl;
  cout << fixed << setprecision(4) << "# sp2OP del_theta:   " << sp2OP_del_theta << endl;
  cout << fixed << setprecision(4) << "# Bin width:         " << binWidth << endl; 
  cout <<                             "# Binning direction: " << binDir << endl; 
 


  ////////////////////////////////////////////////////////////
  // Setting output files
  ////////////////////////////////////////////////////////////

  // For writing coordination number profile of each frame into a file
  const string COORD_NUM_PROFILE = "coord_num_profile_" + TRAJ_FILE + ".dat";
  string rm = "rm -f " + COORD_NUM_PROFILE;
  system(rm.c_str());
  ofstream coordNumProfile(COORD_NUM_PROFILE, ios::app);

  // For writing coordination number profile average (all frames) into a file
  const string COORD_NUM_PROFILE_AVG = "coord_num_profile_avg_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + COORD_NUM_PROFILE_AVG;
  system(rm.c_str());
  ofstream coordNumProfileAvg(COORD_NUM_PROFILE_AVG, ios::app);

  // For writing coordination number profile standard deviation (all frames) into a file
  const string COORD_NUM_PROFILE_STDDEV = "coord_num_profile_stddev_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + COORD_NUM_PROFILE_STDDEV;
  system(rm.c_str());
  ofstream coordNumProfileStdDev(COORD_NUM_PROFILE_STDDEV, ios::app);

  // For writing order parameters profile into a file
  const string OP_PROFILE = "OP_profile_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + OP_PROFILE; 
  system(rm.c_str());
  ofstream OPProfile(OP_PROFILE, ios::app);

  // For writing orderparameter profile average (all frames) into a file
  const string sp2OP_PROFILE_AVG = "sp2OP_profile_avg_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + sp2OP_PROFILE_AVG;
  system(rm.c_str());
  ofstream sp2OPProfileAvg(sp2OP_PROFILE_AVG, ios::app);

  // For writing orderparameter profile average (all frames) into a file
  const string sp3OP_PROFILE_AVG = "sp3OP_profile_avg_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + sp3OP_PROFILE_AVG;
  system(rm.c_str());
  ofstream sp3OPProfileAvg(sp3OP_PROFILE_AVG, ios::app);

  // For writing orderparameter profile standard deviation (all frames) into a file
  const string sp2OP_PROFILE_STDDEV = "sp2OP_profile_stddev_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + sp2OP_PROFILE_STDDEV;
  system(rm.c_str());
  ofstream sp2OPProfileStdDev(sp2OP_PROFILE_STDDEV, ios::app);

  // For writing orderparameter profile standard deviation (all frames) into a file
  const string sp3OP_PROFILE_STDDEV = "sp3OP_profile_stddev_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + sp3OP_PROFILE_STDDEV;
  system(rm.c_str());
  ofstream sp3OPProfileStdDev(sp3OP_PROFILE_STDDEV, ios::app);



  ////////////////////////////////////////////////////////////
  // Reading trajectory file
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
          nbins = traj.boxdims.x / binWidth + 1;
          cout << "# number of bins: " << nbins << endl;
          cout << endl;

          natomsBin.resize(nbins, 0);
          centerBinCoords.resize(nbins, 0.0);
          coordNumBin.resize(nbins, 0.0);
          sp2OPBin.resize(nbins, 0.0);
          sp3OPBin.resize(nbins, 0.0);

          avgCenterBinCoords.resize(nbins, 0.0);
          avgCoordNumBin.resize(nbins, 0.0);
          avgSp2OPBin.resize(nbins, 0.0);
          avgSp3OPBin.resize(nbins, 0.0);

          stdCoordNumBin.resize(nbins, 0.0);
          stdSp2OPBin.resize(nbins, 0.0);
          stdSp3OPBin.resize(nbins, 0.0); 

          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap(); 
        // cout << "# Timestep: " << traj.tstep << endl;
        // cout << "# Lower bounds: " << traj.bounds[0].x << " " << traj.bounds[0].y << " " << traj.bounds[0].z << endl;
        // cout << "# Upper bounds: " << traj.bounds[1].x << " " << traj.bounds[1].y << " " << traj.bounds[1].z << endl;


        // Get neighbor list of each atom in the frame 
        traj.get_neighlist(cutoff);

        // Get coordination number of each atom (rate-determining step)
        // Get order parameters of each atom in the frame
        traj.get_coord_num(cutoff); 
        traj.get_OP(sp3OP_del_theta, sp2OP_del_theta);
        // traj.get_sp3OP(sp3OP_del_theta);
        // traj.get_sp2OP(sp2OP_del_theta);

        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        coordNumBin.assign(coordNumBin.size(), 0);
        sp2OPBin.assign(sp2OPBin.size(), 0);
        sp3OPBin.assign(sp3OPBin.size(), 0);

        // Get and populate order parameters in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].x - traj.bounds[0].x) / binWidth;          
          natomsBin[binIndex]++;

          coordNumBin[binIndex] += traj.coordNums[i];
          sp2OPBin[binIndex] += traj.sp2OP[i];
          sp3OPBin[binIndex] += traj.sp3OP[i];
        }

        // Calculate average coordination number in each bin
        // Calculate average order parameters in each bin
        // Calculate center X coordinate of each bin in the frame
        // Populate avgCoordNumBin and avgCenterBinCoords for all frames
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumBin[i] /= natomsBin[i];
          sp2OPBin[i] /= natomsBin[i];
          sp3OPBin[i] /= natomsBin[i];

          centerBinCoords[i] = traj.bounds[0].x + (i + 0.5) * binWidth;

					if (i == nbins - 1)
					{
					  centerBinCoords[i] = (traj.boxdims.x - i * binWidth) * 0.5 + i * binWidth + traj.bounds[0].x;
					}
          
          // avgCoordNumBin[i] += coordNumBin[i];
          // avgSp2OPBin[i] += sp2OPBin[i];
          // avgSp3OPBin[i] += sp3OPBin[i];
          avgCenterBinCoords[i] += centerBinCoords[i]; 
        }

        // Record coordination number and orderparemeter profiles of the frame for calculation of standard deviation
        coordNumBinAllFrame.push_back(coordNumBin);
        sp2OPBinAllFrame.push_back(sp2OPBin);
        sp3OPBinAllFrame.push_back(sp3OPBin); 
        
        // output order parameters profile into dat file
        // output coordination number profile into dat file
        coordNumProfile << "# timestep: " << traj.tstep << endl; 
        coordNumProfile << "# center_x_coord  avg_coord_num" << endl; 

        OPProfile << "# timestep: " << traj.tstep << endl; 
        OPProfile << "# center_x_coord  avg_sp3OP avg_sp2OP" << endl; 
        
        for (size_t i = 1; i < nbins - 3; i++) // skip the last 2 bins and first bin
        {
          coordNumProfile << fixed << setprecision(4) << centerBinCoords[i] << "  " << coordNumBin[i] << endl;
				  
          OPProfile << fixed << setprecision(4) << centerBinCoords[i] << "  " << sp3OPBin[i] << " " << sp2OPBin[i] << endl;
        }

        nframes++;

        // cout << endl;
      }
      break;

    // For in binning in z direction
    case 2:
      cout << "Binning in z direction" << endl;
      cout << endl;

      while(traj.read_frame()) // Read frames until the end of file
      {
        // Calculate number of bins and initialize vectors
        if (!calledBefore)
        {
          nbins = traj.boxdims.z / binWidth + 1;
          cout << "# number of bins: " << nbins << endl;

          natomsBin.resize(nbins, 0);
          centerBinCoords.resize(nbins, 0.0);
          coordNumBin.resize(nbins, 0.0);
          sp2OPBin.resize(nbins, 0.0);
          sp3OPBin.resize(nbins, 0.0);

          avgCenterBinCoords.resize(nbins, 0.0);
          avgCoordNumBin.resize(nbins, 0.0);
          avgSp2OPBin.resize(nbins, 0.0);
          avgSp3OPBin.resize(nbins, 0.0);

          stdCoordNumBin.resize(nbins, 0.0);
          stdSp2OPBin.resize(nbins, 0.0);
          stdSp3OPBin.resize(nbins, 0.0);

          calledBefore = true; 

        }

        // Wrap coordinates
        traj.pbc_wrap();
        // cout << "# Timestep: " << traj.tstep << endl;
        // cout << "# Lower bounds: " << traj.bounds[0].x << " " << traj.bounds[0].y << " " << traj.bounds[0].z << endl;
        // cout << "# Upper bounds: " << traj.bounds[1].x << " " << traj.bounds[1].y << " " << traj.bounds[1].z << endl;

        // Get neighbor list of each atom in the frame 
        traj.get_neighlist(cutoff);

        // Get coordination number of each atom (rate-determining step)
        // Get order parameters of each atom in the frame
        traj.get_coord_num(cutoff);
        traj.get_OP(sp3OP_del_theta, sp2OP_del_theta); 
        // traj.get_sp3OP(sp3OP_del_theta);
        // traj.get_sp2OP(sp2OP_del_theta);
      
        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        coordNumBin.assign(coordNumBin.size(), 0);
        sp2OPBin.assign(sp2OPBin.size(), 0);
        sp3OPBin.assign(sp3OPBin.size(), 0);

        // Get and populate order parameters in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].z - traj.bounds[0].z) / binWidth;
          natomsBin[binIndex]++;

          coordNumBin[binIndex] += traj.coordNums[i];
          sp2OPBin[binIndex] += traj.sp2OP[i];
          sp3OPBin[binIndex] += traj.sp3OP[i];
        }

        // Calculate average coordination number in each bin
        // Calculate average order parameters in each bin
        // Calculate center Z coordinate of each bin in the frame
        // Populate avgCoordNumBin and avgCenterBinCoords for all frames
        for (size_t i = 0; i < nbins; i++)
        {
          coordNumBin[i] /= natomsBin[i];
          sp2OPBin[i] /= natomsBin[i];
          sp3OPBin[i] /= natomsBin[i];

          centerBinCoords[i] = traj.bounds[0].z + (i + 0.5) * binWidth;

					if (i == nbins - 1)
					{
					  centerBinCoords[i] = (traj.boxdims.z - i * binWidth) * 0.5 + i * binWidth + traj.bounds[0].z;
					}
          
          // avgCoordNumBin[i] += coordNumBin[i];
          // avgSp2OPBin[i] += sp2OPBin[i];
          // avgSp3OPBin[i] += sp3OPBin[i];
          avgCenterBinCoords[i] += centerBinCoords[i]; 
        }

        // Record coordination number and orderparemeter profiles of the frame for calculation of standard deviation
        coordNumBinAllFrame.push_back(coordNumBin);
        sp2OPBinAllFrame.push_back(sp2OPBin);
        sp3OPBinAllFrame.push_back(sp3OPBin); 
        
        // output coordination number profile into dat file
        // output order parameters profile into dat file
        coordNumProfile << "# timestep: " << traj.tstep << endl; 
        coordNumProfile << "# center_z_coord  avg_coord_num" << endl; 

        OPProfile << "# timestep: " << traj.tstep << endl; 
        OPProfile << "# center_z_coord  avg_sp3OP avg_sp2OP" << endl; 

        for (size_t i = 1; i < nbins - 3; i++)  // skip the last 2 bins and first bin
        {
          coordNumProfile << fixed << setprecision(4) << centerBinCoords[i] << "  " << coordNumBin[i] << endl;

				  OPProfile << fixed << setprecision(4) << centerBinCoords[i] << "  " << sp3OPBin[i] << " " << sp2OPBin[i] << endl;				
        }

        nframes++;

        // cout << endl;
      }
      break;
  }

  // Calculate average of the quantities of each bin among all frames
  // Calculate the average center coordinate of each bin among all frames
  for (size_t i = 0; i < nbins; i++)
  {
    // avgCoordNumBin[i] /= nframes;
    // avgSp2OPBin[i] /= nframes;
    // avgSp3OPBin[i] /= nframes;
    avgCenterBinCoords[i] /= nframes; 
  }

  // Calculate average and standard deviation of the quantities of each bin among all frames
  // Note: dimension of the "all-frame" vectors: (nrow = nframes, ncol = nbins) 
  vector<double> BinCoordNumOfAllFrame; // column vector of the all-frame vector
  vector<double> BinSp2OpOfAllFrame;  // column vector of the all-frame vector
  vector<double> BinSp3OpOfAllFrame;  // column vector of the all-frame vector

  BinCoordNumOfAllFrame.resize(nframes, 0.0); 
  BinSp2OpOfAllFrame.resize(nframes, 0.0); 
  BinSp3OpOfAllFrame.resize(nframes, 0.0); 

  for (size_t j = 0; j < nbins; j++)
  {
    for (size_t i = 0; i < nframes; i++)
    {
      BinCoordNumOfAllFrame[i] = coordNumBinAllFrame[i][j];
      BinSp2OpOfAllFrame[i]    = sp2OPBinAllFrame[i][j];
      BinSp3OpOfAllFrame[i]    = sp3OPBinAllFrame[i][j];
    }

    avgCoordNumBin[j] = mean(BinCoordNumOfAllFrame);
    avgSp2OPBin[j]    = mean(BinSp2OpOfAllFrame);
    avgSp3OPBin[j]    = mean(BinSp3OpOfAllFrame);

    stdCoordNumBin[j] = stddev(BinCoordNumOfAllFrame);
    stdSp2OPBin[j]    = stddev(BinSp2OpOfAllFrame);
    stdSp3OPBin[j]    = stddev(BinSp3OpOfAllFrame);
  }


  // Output average coordination number profile
  // Output stddev coordination number profile 
  // Output average orderparmeter profile
  // Output stddev orderparmeter profile 
  coordNumProfileAvg << "# number of timesteps averaged: " << nframes << endl; 
  coordNumProfileAvg << "# avg_center_coord  avg_coord_num_all_frames" << endl;

  coordNumProfileStdDev << "# number of timesteps averaged: " << nframes << endl; 
  coordNumProfileStdDev << "# avg_center_coord  stddev_coord_num_all_frames" << endl;

  sp2OPProfileAvg << "# number of timesteps averaged: " << nframes << endl; 
  sp2OPProfileAvg << "# avg_center_coord  avg_sp2OP_all_frames" << endl;

  sp3OPProfileAvg << "# number of timesteps averaged: " << nframes << endl; 
  sp3OPProfileAvg << "# avg_center_coord  avg_sp3OP_all_frames" << endl;

  sp2OPProfileStdDev << "# number of timesteps averaged: " << nframes << endl; 
  sp2OPProfileStdDev << "# avg_center_coord  stddev_sp2OP_all_frames" << endl; 

  sp3OPProfileStdDev << "# number of timesteps averaged: " << nframes << endl; 
  sp3OPProfileStdDev << "# avg_center_coord  stddev_sp3OP_all_frames" << endl;   

  for (size_t i = 0; i < nbins - 2; i++) // skip the last 2 bins (adjust as needed)
  {
    coordNumProfileAvg << fixed << setprecision(4) << avgCenterBinCoords[i] << "  " << avgCoordNumBin[i] << endl;

    coordNumProfileStdDev << fixed << setprecision(4) << avgCenterBinCoords[i] << "  " << stdCoordNumBin[i] << endl; 
   
    sp2OPProfileAvg << fixed << setprecision(4) << avgCenterBinCoords[i] << "  " << avgSp2OPBin[i] << endl;

    sp3OPProfileAvg << fixed << setprecision(4) << avgCenterBinCoords[i] << "  " << avgSp3OPBin[i] << endl;

    sp2OPProfileStdDev << fixed << setprecision(4) << avgCenterBinCoords[i] << "  " << stdSp2OPBin[i] << endl; 

    sp3OPProfileStdDev << fixed << setprecision(4) << avgCenterBinCoords[i] << "  " << stdSp3OPBin[i] << endl; 
  }
 

  // End frame processing
  cout << endl;
  cout << "# " << nframes  << " frames are read and processed" << endl;
  cout << endl;


  // Close files and confirm writing
  coordNumProfile.close();
  cout << "# Coordination number profile written to " << COORD_NUM_PROFILE << endl;
  // cout << endl;

  coordNumProfileAvg.close();
  cout << "# Average coordination number profile written to " << COORD_NUM_PROFILE_AVG << endl;
  // cout << endl;

  coordNumProfileStdDev.close();
  cout << "# Standard deviation of coordination number profile is written to " << COORD_NUM_PROFILE_STDDEV << endl;
  // cout << endl;


  OPProfile.close();
  cout << "# Order parameter profiles are written to " << OP_PROFILE << endl;
  // cout << endl;
  
  sp2OPProfileAvg.close();
  cout << "# Average sp2 order parameter profile is written to " << sp2OP_PROFILE_AVG << endl;
  // cout << endl;

  sp3OPProfileAvg.close();
  cout << "# Average sp3 order parameter profile is written to " << sp3OP_PROFILE_AVG << endl;
  // cout << endl;

  sp2OPProfileStdDev.close();
  cout << "# Standard deviation of sp2 order parameter profile is written to " << sp2OP_PROFILE_STDDEV << endl;


  sp3OPProfileStdDev.close();
  cout << "# Standard deviation of sp3 order parameter profile is written to " << sp3OP_PROFILE_STDDEV << endl;
 

  return 0; 
}
