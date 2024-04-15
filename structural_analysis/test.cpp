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
  // Exit execution if input argument is incorrect
  if (argc != 4)
    {
      cerr << "Usage: ./<this_file> <traj_file> <neighlist_cutoff> <sp3OP_delta_theta>\n";   // cerr will not be flushed.
	    return 1;
    }
 
  const string TRAJ_FILE = argv[1];  // Only 4 atoms in the test traj file
  Trajectory traj(TRAJ_FILE); // Create the traj object

  // Cutoff distance for building neighbor list and calculating sp3OP
  double cutoff = stod(argv[2]);
  // Tolerance angle theta for sp3 order parameter
  double del_theta = stod(argv[3]);

  traj.read_frame();

  traj.pbc_wrap();

  traj.get_neighlist(cutoff);

  traj.get_sp3OP(del_theta);
  traj.get_sp2OP(del_theta);

  // // for debugging neighlist
  // cout << endl;
  // for (size_t i = 0; i < traj.natoms; i++)
  // {
  //   cout << "Neighbor list of atom " << i << ": "; 
  //   for (size_t j = 0; j < traj.neighlist[i].size(); j++)
  //   {
  //     cout << traj.neighlist[i][j] << " "; 
  //   }
  //   cout << endl;
  //   cout << fixed << setprecision(4) << "sp3OP of atom: " << i << ": " << traj.sp3OP[i] << endl; 
  //   cout << fixed << setprecision(4) << "sp2OP of atom: " << i << ": " << traj.sp2OP[i] << endl;
  //   cout << endl;
  // }  

  return 0;
}