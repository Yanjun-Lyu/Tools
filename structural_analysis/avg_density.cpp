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
  if (argc != 8)
    {
      cerr << "Usage: ./<this_file> <traj_file> xlo xhi ylo yhi zlo zhi \n";   // cerr will not be flushed.
	    return 1;
    }

  const string TRAJ_FILE = argv[1];
  Trajectory traj(TRAJ_FILE);

  ////////////////////////////////////////////////////////////
  // Reading trajectory file
  ////////////////////////////////////////////////////////////

  int nframes = 1;
  bool calledBefore = false;

  // For user-specified bounds vector
  double xlo = stod(argv[2]);
  double xhi = stod(argv[3]);
  double ylo = stod(argv[4]);
  double yhi = stod(argv[5]);
  double zlo = stod(argv[6]);
  double zhi = stod(argv[7]);

  vector<xyz> us_sp_bounds;
  us_sp_bounds.push_back({xlo, ylo, zlo});
  us_sp_bounds.push_back({xhi, yhi, zhi}); 

  cout << fixed << setprecision(4) << "# Calculate density in (xlo xhi ylo yhi zlo zhi): " << xlo << " " << xhi << " " << ylo << " " << yhi << " " << zlo << " " << zhi << endl;

  bool isInX;
  bool isInY;
  bool isInZ;
  bool isInXYZ;

  traj.read_frame(); // Read a single frame!

  // check if us_sp_bound is within box dimensions
  isInX = xlo >= traj.bounds[0].x && xhi <= traj.bounds[1].x;
  isInY = ylo >= traj.bounds[0].y && yhi <= traj.bounds[1].y;
  isInZ = zlo >= traj.bounds[0].z && zhi <= traj.bounds[1].z;
  isInXYZ = isInX && isInY && isInZ;

  if (!isInXYZ)
  {
    cerr << "# User-specified bounds out of range.";
    return 1;
  }
      
  // Wrap coordinates
  traj.pbc_wrap(); 

  double density = traj.get_mass_den(us_sp_bounds);

  cout << fixed << setprecision(4) << "# Density: " << density << " gcc" << endl;

  return 0;
}