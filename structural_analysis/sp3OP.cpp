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
  if (argc != 6)
    {
      cerr << "Usage: ./<this_file> <traj_file> <neighlist_cutoff> <sp3OP_delta_theta> <profile_bin_width> <binning_direction (1 for x, 2 for z)>\n";   // cerr will not be flushed.
	    return 1;
    }


  ////////////////////////////////////////////////////////////
  // Setting parameters
  ////////////////////////////////////////////////////////////

  const string TRAJ_FILE = argv[1];  // Only 4 atoms in the test traj file
  Trajectory traj(TRAJ_FILE); // Create the traj object

  // Cutoff distance for building neighbor list and calculating sp3OP
  double cutoff = stod(argv[2]);
  // Tolerance angle theta for sp3 order parameter
  double sp3OP_del_theta = stod(argv[3]);

  // Bins for sp3 orderparameter profile
  int nbins;                           // Number of bins
  int binIndex;                        // Bin index of an atom of interest
  double binWidth = stod(argv[4]);     // Bin width for sp3OP profile (Ang)

  unsigned int binDir = stoi(argv[5]); // Binning direction (x, y or z)

  vector<int> natomsBin;        // Number of atoms in the bin of interest
  vector<double> sp3OPBin;   // total coordination number of atom in the bin


  ////////////////////////////////////////////////////////////
  // Setting output files
  ////////////////////////////////////////////////////////////

  // For writing sp3 order parameter profile into a file
  const string SP3OP_PROFILE = "sp3OP_profile_" + TRAJ_FILE + ".dat";
  string rm = "rm -f " + SP3OP_PROFILE; 
  system(rm.c_str());
  ofstream sp3Profile(SP3OP_PROFILE, ios::app);

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
          sp3OPBin.resize(nbins, 0.0);
          
          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap(); 

        // Get neighbor list of each atom in the frame 
        traj.get_neighlist(cutoff);

        // Get sp3 order parameter of each atom in the frame
        traj.get_sp3OP(sp3OP_del_theta);
     
        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        sp3OPBin.assign(sp3OPBin.size(), 0);

        // Get and populate sp3 order parameter in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].x - traj.bounds[0].x) / binWidth;
          natomsBin[binIndex]++;
          sp3OPBin[binIndex] += traj.sp3OP[i];
        }

        // Calculate average sp3 order parameter in each bin
        for (size_t i = 0; i < nbins; i++)
        {
          sp3OPBin[i] /= natomsBin[i];
        }
        
        // output sp3 order parameter profile into dat file
        sp3Profile << "# timestep: " << traj.tstep << endl; 
        sp3Profile << "# center_x_coord  avg_sp3OP" << endl; 
        for (size_t i = 0; i < nbins; i++)
        {
          sp3Profile << fixed << setprecision(4) << traj.bounds[0].x + (i - 0.5) * binWidth << "  " << sp3OPBin[i] << endl;
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
          nbins = traj.boxdims.z / binWidth;
          cout << "# number of bins: " << nbins << endl;

          natomsBin.resize(nbins,0);
          sp3OPBin.resize(nbins, 0.0);
          
          calledBefore = true;
        }

        // Wrap coordinates
        traj.pbc_wrap();

        // Get neighbor list of each atom in the frame 
        traj.get_neighlist(cutoff);

        // Get sp3 order parameter of each atom in the frame
        traj.get_sp3OP(sp3OP_del_theta);
      
        // Initialize vars
        natomsBin.assign(natomsBin.size(), 0);
        sp3OPBin.assign(sp3OPBin.size(), 0);

        // Get and populate sp3 order parameter in each bin
        for (size_t i = 0; i < traj.natoms; i++)
        {
          binIndex = (traj.coords[i].z - traj.bounds[0].z) / binWidth;
          natomsBin[binIndex]++;
          sp3OPBin[binIndex] += traj.sp3OP[i];
        }

        // Calculate average sp3 order parameter in each bin
        for (size_t i = 0; i < nbins; i++)
        {
          sp3OPBin[i] /= natomsBin[i];
        }
        
        // output sp3 order parameter profile into dat file
        sp3Profile << "# timestep: " << traj.tstep << endl; 
        sp3Profile << "# center_z_coord  avg_sp3OP" << endl; 
        for (size_t i = 0; i < nbins; i++)
        {
          sp3Profile << fixed << setprecision(4) << traj.bounds[0].z + (i - 0.5) * binWidth << "  " << sp3OPBin[i] << endl;
        }

        nframes++;

        cout << endl;
      }
      break;
  }
  
  cout << endl;
  cout << "# " << nframes - 1  << " frames are read and processed" << endl;
  cout << endl;

  // Close file and confirm writing.
  sp3Profile.close();
  cout << "# sp3 order parameter profile written to " << SP3OP_PROFILE << endl;
  
  return 0;
}






  // sp2OP
  // double sp2OP_del_theta = stod(argv[4]);
  // traj.get_sp2OP(sp2OP_del_theta);
  
  // for debugging neighlist
  // cout << endl;
  // for (size_t i = 0; i < traj.natoms; i++)
  // {
  //   cout << "Neighbor list of atom: " << i << ": "; 
  //   for (size_t j = 0; j < traj.neighlist[i].size(); j++)
  //   {
  //     cout << traj.neighlist[i][j] << " "; 
  //   }
  //   cout << endl;
  // }  


  // cout << fixed << setprecision(4) << "The sp3OP of atom 0: " << traj.sp3OP[0] << endl;
 

  // // test get_azimuth (traj_90degree.lammpstrj)
  // xyz ij, ik, il, im; // vectors for testing vector and angle calculations
  // ij = traj.get_vec(0, 1);
  // cout << "01: " << ij.x << " " << ij.y << " " << ij.z << endl;
  // ik = traj.get_vec(0, 2);
  // cout << "02: " << ik.x << " " << ik.y << " " << ik.z << endl;
  // il = traj.get_vec(0, 3);
  // cout << "03: " << il.x << " " << il.y << " " << il.z << endl;
  // im = traj.get_vec(0, 4);
  // cout << "04: " << im.x << " " << im.y << " " << im.z << endl;
 
  // double angle_ij_ik = traj.get_angle(ij, ik);
  // double angle_ij_il = traj.get_angle(ij, il);
  // double angle_ij_im = traj.get_angle(ij, im);
  // double angle_ik_il = traj.get_angle(ik, il);
  // double angle_ik_im = traj.get_angle(ik, im);
  // double angle_il_im = traj.get_angle(il, im);

  // cout << "# Angle      0102: "  << angle_ij_ik     << endl;
  // cout << "# Angle      0103: "  << angle_ij_il     << endl;
  // cout << "# Angle      0104: "  << angle_ij_im     << endl;
  // cout << "# Angle      0203: "  << angle_ik_il     << endl;
  // cout << "# Angle      0204: "  << angle_ik_im     << endl;
  // cout << "# Angle      0304: "  << angle_il_im     << endl;

  // double azimuth_ijik_il = traj.get_azimuth(ij, ik, il);
  // double azimuth_ijim_il = traj.get_azimuth(ij, im, il);
  // double azimuth_ikim_il = traj.get_azimuth(ik, im, il);
  
  // double azimuth_ijik_im = traj.get_azimuth(ij, ik, im);
  // double azimuth_ijil_im = traj.get_azimuth(ij, il, im); 
  // double azimuth_ikil_im = traj.get_azimuth(ik, il, im);

  // double azimuth_ijil_ik = traj.get_azimuth(ij, il, ik);
  // double azimuth_ijim_ik = traj.get_azimuth(ij, im, ik);
  // double azimuth_ilim_ik = traj.get_azimuth(il, im, ik); 

  // double azimuth_ikil_ij = traj.get_azimuth(ik, il, ij);
  // double azimuth_ikim_ij = traj.get_azimuth(ik, im, ij);
  // double azimuth_ilim_ij = traj.get_azimuth(il, im, ij);

  // cout << "# Azimuth 0102_03: "  << azimuth_ijik_il << endl; 
  // cout << "# Azimuth 0104_03: "  << azimuth_ijim_il << endl; 
  // cout << "# Azimuth 0204_03: "  << azimuth_ikim_il << endl;

  // cout << "# Azimuth 0102_04: "  << azimuth_ijik_im << endl; 
  // cout << "# Azimuth 0103_04: "  << azimuth_ijil_im << endl; 
  // cout << "# Azimuth 0203_04: "  << azimuth_ikil_im << endl; 

  // cout << "# Azimuth 0103_02: "  << azimuth_ijil_ik << endl; 
  // cout << "# Azimuth 0104_02: "  << azimuth_ijim_ik << endl; 
  // cout << "# Azimuth 0304_02: "  << azimuth_ilim_ik << endl;  

  // cout << "# Azimuth 0203_01: "  << azimuth_ikil_ij << endl; 
  // cout << "# Azimuth 0204_01: "  << azimuth_ikim_ij << endl; 
  // cout << "# Azimuth 0304_01: "  << azimuth_ilim_ij << endl; 


 // cout << fixed << setprecision(4) << "i: " << traj.coords[0].x << " " << traj.coords[0].y << " " << traj.coords[0].z << endl;


  // cout << fixed << setprecision(4) << "ij: " << ij.x << " " << ij.y << " " << ij.z << endl; 
  // cout << fixed << setprecision(4) << "ik: " << ik.x << " " << ik.y << " " << ik.z << endl; 
  // cout << fixed << setprecision(4) << "il: " << il.x << " " << il.y << " " << il.z << endl; 

  // cout << fixed << setprecision(4) << "rij: " << l2norm(ij) << endl; 
  // cout << fixed << setprecision(4) << "rik: " << l2norm(ik) << endl; 
  // cout << fixed << setprecision(4) << "ril: " << l2norm(il) << endl; 

  // cout << fixed << setprecision(4) << "dot_ij_ik: " << dot(ij, ik) << endl;
  // cout << fixed << setprecision(4) << "dot_ij_il: " << dot(ij, il) << endl;
  // cout << fixed << setprecision(4) << "dot_ik_il: " << dot(ik, il) << endl;