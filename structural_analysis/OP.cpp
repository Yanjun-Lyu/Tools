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
      cerr << "Usage: ./<this_file> <traj_file> <neighlist_cutoff> <sp3OP_delta_theta> <sp2OP_delta_theta> <profile_bin_width> <binning_direction (1 for x, 2 for z)>\n";   // cerr will not be flushed.
	    return 1;
    }


  ////////////////////////////////////////////////////////////
  // Setting parameters
  ////////////////////////////////////////////////////////////

  const string TRAJ_FILE = argv[1];  // Only 4 atoms in the test traj file
  Trajectory traj(TRAJ_FILE); // Create the traj object

  // Cutoff distance for building neighbor list and calculating sp2OP
  double cutoff = stod(argv[2]);
  // Tolerance angle theta for order parameters
  double sp3OP_del_theta = stod(argv[3]);
  double sp2OP_del_theta = stod(argv[4]);

  // Bins for orderparameter profile
  int nbins;                           // Number of bins
  int binIndex;                        // Bin index of an atom of interest
  double binWidth = stod(argv[5]);     // Bin width for sp2OP profile (Ang)

  unsigned int binDir = stoi(argv[6]); // Binning direction (x, y or z)

  vector<int> natomsBin;        // Number of atoms in the bin of interest
  vector<double> sp3OPBin;      // total sp3OP of atom in the bin (populated)
  vector<double> sp2OPBin;      // total sp2OP of atom in the bin (populated)

  //for debugging
  //vector<vector<int>> atomsInBin;  // Which atoms are in each bin

  vector<double> centerXCoords;  // Average x coordinate of each bin index in one frame
  vector<double> avgCenterXCoords;    // Average x coordinate of each bin index in all frames
  vector<double> sp2OPBinAllFrame;   // total sp2OP of atoms in the bin in all frames
  vector<double> sp3OPBinAllFrame;   // total sp3OP of atoms in the bin in all frames


  cout << "# Input parameters:" << endl;
  cout << fixed << setprecision(4) << "# neighlist cutoff:  " << cutoff << endl; 
  cout << fixed << setprecision(4) << "# sp3OP del_theta:   " << sp3OP_del_theta << endl;
  cout << fixed << setprecision(4) << "# sp2OP del_theta:   " << sp2OP_del_theta << endl;
  cout << fixed << setprecision(4) << "# Bin width:         " << binWidth << endl; 
  cout <<                             "# Binning direction: " << binDir << endl; 
   
  ////////////////////////////////////////////////////////////
  // Setting output files
  ////////////////////////////////////////////////////////////

  // For writing order parameters profile into a file
  const string OP_PROFILE = "OP_profile_" + TRAJ_FILE + ".dat";
  string rm = "rm -f " + OP_PROFILE; 
  system(rm.c_str());
  ofstream OPProfile(OP_PROFILE, ios::app);

  // For writing coordination number profile average into a file
  const string OP_AVG = "OP_avg_" + TRAJ_FILE + ".dat";
  rm = "rm -f " + OP_AVG;
  system(rm.c_str());
  ofstream OPProfileAvg(OP_AVG, ios::app);


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

          natomsBin.resize(nbins, 0);
          sp2OPBin.resize(nbins, 0.0);
          sp3OPBin.resize(nbins, 0.0);
          centerXCoords.resize(nbins, 0.0);
          avgCenterXCoords.resize(nbins, 0.0);
          sp2OPBinAllFrame.resize(nbins, 0.0);
          sp3OPBinAllFrame.resize(nbins, 0.0);

          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap(); 
        // cout << "# Timestep: " << traj.tstep << endl;
        // cout << "# Lower bounds: " << traj.bounds[0].x << " " << traj.bounds[0].y << " " << traj.bounds[0].z << endl;
        // cout << "# Upper bounds: " << traj.bounds[1].x << " " << traj.bounds[1].y << " " << traj.bounds[1].z << endl;


        // Get neighbor list of each atom in the frame 
        traj.get_neighlist(cutoff);

        // Get order parameters of each atom in the frame
        traj.get_sp3OP(sp3OP_del_theta);
        traj.get_sp2OP(sp2OP_del_theta);
     
        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        sp2OPBin.assign(sp2OPBin.size(), 0);
        sp3OPBin.assign(sp3OPBin.size(), 0);

        // Get and populate order parameters in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].x - traj.bounds[0].x) / binWidth;          
          natomsBin[binIndex]++;
          sp2OPBin[binIndex] += traj.sp2OP[i];
          sp3OPBin[binIndex] += traj.sp3OP[i];
        }


        // Calculate average order parameters in each bin
        // Calculate center X coordinate number of each bin in the frame
        // Populate coordNumBinAllFrame and avgCenterXCoords for all frames
        for (size_t i = 0; i < nbins; i++)
        {
          sp2OPBin[i] /= natomsBin[i];
          sp3OPBin[i] /= natomsBin[i];

          centerXCoords[i] = traj.bounds[0].x + (i + 0.5) * binWidth;
					if (i == nbins - 1)
					{
					  centerXCoords[i] = (traj.boxdims.x - i * binWidth) * 0.5 + i * binWidth + traj.bounds[0].x;
					}
          
          sp2OPBinAllFrame[i] += sp2OPBin[i];
          sp3OPBinAllFrame[i] += sp3OPBin[i];
          avgCenterXCoords[i] += centerXCoords[i]; 

        }
        
        // output order parameters profile into dat file
        OPProfile << "# timestep: " << traj.tstep << endl; 
        OPProfile << "# center_x_coord  avg_sp3OP avg_sp2OP" << endl; 
        for (size_t i = 0; i < nbins; i++) 
        {
				  OPProfile << fixed << setprecision(4) << centerXCoords[i] << "  " << sp3OPBin[i] << " " << sp2OPBin[i] << endl;
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

        // Calculate number of bins and initialize vectors
        if (!calledBefore)
        {
          nbins = traj.boxdims.z / binWidth + 1;
          cout << "# number of bins: " << nbins << endl;

          natomsBin.resize(nbins,0);
          sp2OPBin.resize(nbins, 0.0);
          sp3OPBin.resize(nbins, 0.0);
          centerXCoords.resize(nbins, 0.0);
          avgCenterXCoords.resize(nbins, 0.0);
          sp2OPBinAllFrame.resize(nbins, 0.0);
          sp3OPBinAllFrame.resize(nbins, 0.0);
          
          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap();
        // cout << "# Timestep: " << traj.tstep << endl;
        // cout << "# Lower bounds: " << traj.bounds[0].x << " " << traj.bounds[0].y << " " << traj.bounds[0].z << endl;
        // cout << "# Upper bounds: " << traj.bounds[1].x << " " << traj.bounds[1].y << " " << traj.bounds[1].z << endl;

        // Get neighbor list of each atom in the frame 
        traj.get_neighlist(cutoff);

        // Get order parameters of each atom in the frame
        traj.get_sp3OP(sp3OP_del_theta);
        traj.get_sp2OP(sp2OP_del_theta);
      
        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        sp2OPBin.assign(sp2OPBin.size(), 0);
        sp3OPBin.assign(sp3OPBin.size(), 0);

        // Get and populate order parameters in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].z - traj.bounds[0].z) / binWidth;
          natomsBin[binIndex]++;
          sp2OPBin[binIndex] += traj.sp2OP[i];
          sp3OPBin[binIndex] += traj.sp3OP[i];
        }

        // Calculate average order parameters in each bin
        // Calculate center Z coordinate number of each bin in the frame
        // Populate coordNumBinAllFrame and avgCenterXCoords for all frames
        for (size_t i = 0; i < nbins; i++)
        {
          sp2OPBin[i] /= natomsBin[i];
          sp3OPBin[i] /= natomsBin[i];

          centerXCoords[i] = traj.bounds[0].z + (i + 0.5) * binWidth;
					if (i == nbins - 1)
					{
					  centerXCoords[i] = (traj.boxdims.z - i * binWidth) * 0.5 + i * binWidth + traj.bounds[0].z;
					}
          
          sp2OPBinAllFrame[i] += sp2OPBin[i];
          sp3OPBinAllFrame[i] += sp3OPBin[i];
          avgCenterXCoords[i] += centerXCoords[i]; 

        }
        
        // output order parameters profile into dat file
        OPProfile << "# timestep: " << traj.tstep << endl; 
        OPProfile << "# center_z_coord  avg_sp3OP avg_sp2OP" << endl; 
        for (size_t i = 0; i < nbins; i++)
        {
				  OPProfile << fixed << setprecision(4) << centerXCoords[i] << "  " << sp3OPBin[i] << " " << sp2OPBin[i] << endl;				
        }

        nframes++;

        cout << endl;
      }
      break;
  }

  // Calculate average coordination number and center x coordinates among all frames
  for (size_t i = 0; i < nbins; i++)
  {
    sp2OPBinAllFrame[i] /= nframes;
    sp3OPBinAllFrame[i] /= nframes;
    avgCenterXCoords[i] /= nframes; 
  }

  // Output average coordination number profile into dat file
  OPProfileAvg << "# number of timesteps averaged: " << nframes << endl; 
  OPProfileAvg << "# avg_center_coord  avg_sp3OP_all_frames avg_sp2OP_all_frames" << endl; 
  for (size_t i = 0; i < nbins; i++)
  {
    OPProfileAvg << fixed << setprecision(4) << avgCenterXCoords[i] << "  " << sp3OPBinAllFrame[i] << " " << sp2OPBinAllFrame[i] << endl;
  }
 
  cout << endl;
  cout << "# " << nframes  << " frames are read and processed" << endl;
  cout << endl;

  // Close file and confirm writing.
  OPProfile.close();
  cout << "# Order parameter profiles are written to " << OP_PROFILE << endl;
  
  OPProfileAvg.close();
  cout << "# Average order parameter profile written to " << OP_AVG << endl;
  return 0;
}
