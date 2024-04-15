#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

#include "helpers.hpp"
#include "trajectory.hpp"

using namespace std;


int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "Usage: ./this_file traj_file\n";
	return 1;
    }
    
    Trajectory traj;
    
    xyz A;
    tokenize(get_next_line());
    
    
    
    
    A.x = stod(argv[1]);
    A.y = stod(argv[2]);
    A.z = stod(argv[3]);

    xyz B;
    B.x = stod(argv[4]);
    B.y = stod(argv[5]);
    B.z = stod(argv[6]);    
    
    traj.frame.push_back({A});
    traj.frame.push_back({B}); 
   
    traj.boxdim.x = 5;
    traj.boxdim.y = 5;
    traj.boxdim.z = 5;
             
    double dist = traj.get_pbc_dist(traj.frame[0], traj.frame[1], traj.boxdim);
         
    traj.frame = traj.pbc_wrap(traj.frame, traj.boxdim);
    
    cout << "Number of atoms: " << traj.frame.size() << endl;
    cout << "Wrapped Coordinates of the first atom: " << traj.frame[0].x << ", " << traj.frame[0].y << ", " << traj.frame[0].z << endl; 
    cout << "Wrapped Coordinates of the second atom: " << traj.frame[1].x << ", " << traj.frame[1].y << ", " << traj.frame[1].z << endl;
    cout << "Distance between the first and second atoms: " << dist << endl;
    
    return 0;
}
