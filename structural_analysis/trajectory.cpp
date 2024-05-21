#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

#include "trajectory.hpp"
#include "helpers.hpp"

using namespace std;

Trajectory::Trajectory(string trajFile):bounds(2)
{
    cout << "# Read from trajectory file: " << trajFile << endl;
    trajstream.open(trajFile);
}


Trajectory::~Trajectory()
{
    trajstream.close();
}


// Read frame in lammpstrj file
bool Trajectory::read_frame()
{    
    string  line;
    static bool called_before = false;
    
    // Ensure we haven't reached the end of the trajectory
    if (!get_next_line(trajstream, line)) // ITEM: TIMESTEP (if cannot get it) 
    {
        cout << "# End of trajectory file detected" << endl;
            return false;
    }

    // Read current timestep
    tstep = stoi(get_next_line(trajstream));  
    cout << "# Reading timestep: " << tstep << endl;

    get_next_line(trajstream);  // ITEM: NUMBER OF ATOMS
    
    if (!called_before)
    {
        // Read the number of atoms   
        natoms = stoi(get_next_line(trajstream));
        cout << "# Number of atoms: " << natoms << endl;
    }
    else
    {
        get_next_line(trajstream);  // natoms 
    }

    // Read the box dimensions; assume they start from 0,0,0 (for NVT)
    get_next_line(trajstream);  // ITEM: BOX BOUNDS pp pp pp
    
    line = get_next_line(trajstream); 
    bounds[0].x = stod(get_token(line, 0));
    bounds[1].x = stod(get_token(line, 1));

    line = get_next_line(trajstream); 
    bounds[0].y = stod(get_token(line, 0));
    bounds[1].y = stod(get_token(line, 1));
    
    line = get_next_line(trajstream); 
    bounds[0].z = stod(get_token(line, 0));
    bounds[1].z = stod(get_token(line, 1));
    
    boxdims.x = bounds[1].x - bounds[0].x;
    boxdims.y = bounds[1].y - bounds[0].y; 
    boxdims.z = bounds[1].z - bounds[0].z; 

    // cout << "# Box dimensions: " << boxdims.x << " " << boxdims.y << " " << boxdims.z << endl;

    // // Calculate number density
    // numDen = natoms/(boxdims.x)/(boxdims.y)/(boxdims.z);
    
    // Read the configuration coordinates
    get_next_line(trajstream);  // ITEM: ATOMS id type element xu yu zu
    
    coords.clear();  // Size: 0 (no element);   Capacity: dynamic.
    xyz     coordinate;
    
    for (size_t i = 0; i < natoms; i++)
    {
        line = get_next_line(trajstream);
        
        coordinate.x = stod(get_token(line,3));
        coordinate.y = stod(get_token(line,4));
        coordinate.z = stod(get_token(line,5));
	
        coords.push_back(coordinate);
    }
    
    if (!called_before)
        called_before = true;
    
    
    return true;
}


// Wrap atom coordinates according to periodic boundary condition
void Trajectory::pbc_wrap()
{
    for (size_t i = 0; i < natoms; i++)  // size_t is equivalent to unsigned int for presenting the size of objects.
    {
        coords[i].x -= boxdims.x * floor((coords[i].x - bounds[0].x) / boxdims.x);
        coords[i].y -= boxdims.y * floor((coords[i].y - bounds[0].y) / boxdims.y);
        coords[i].z -= boxdims.z * floor((coords[i].z - bounds[0].z) / boxdims.z);
    }
}


// Get distance between two atoms, applied periodic boundary condition
double Trajectory::get_dist(int i, int j)
{
    double dx = coords[j].x - coords[i].x;
    double dy = coords[j].y - coords[i].y;
    double dz = coords[j].z - coords[i].z;
    
    dx -= boxdims.x * round(dx / boxdims.x);
    dy -= boxdims.y * round(dy / boxdims.y);
    dz -= boxdims.z * round(dz / boxdims.z);

    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}



// Return density within a user-specified square area in g/cm^3
double Trajectory::get_mass_den(vector<xyz> usr_sp_bounds)
{
    const double MASS_C_ATOM = 1.9944733e-23;         // A carbon atom mass in g
    const double A2CM = 1e-8;                         // Angstrom to centimeter

    // Caculate the user-specified area
    double dx = usr_sp_bounds[1].x - usr_sp_bounds[0].x;
    double dy = usr_sp_bounds[1].y - usr_sp_bounds[0].y;
    double dz = usr_sp_bounds[1].z - usr_sp_bounds[0].z;
    double area = dx * A2CM * dy * A2CM * dz * A2CM;  // area in cm^3
    
    // Count number of atoms in user-specified area
    bool inX, inY, inZ;
    int count = 0;
    for (size_t i = 0; i < natoms; i++)
    {
        inX = coords[i].x >= usr_sp_bounds[0].x && coords[i].x <= usr_sp_bounds[1].x;
        inY = coords[i].y >= usr_sp_bounds[0].y && coords[i].y <= usr_sp_bounds[1].y; 
        inZ = coords[i].z >= usr_sp_bounds[0].z && coords[i].z <= usr_sp_bounds[1].z; 
        if (inX && inY && inZ)
        {
            count++;
        }
    }
    cout << "# Number of atoms in user-specified bounds: " << count << endl;
    // Calculate density in g/cm^3
    double massDen = count * MASS_C_ATOM / area;
    return massDen;
}



// Calculate coordination number of each atom
void Trajectory::get_coord_num(double cutoff)
{
    coordNums.clear();
    coordNums.resize(natoms, 0);
      
    double dist;
    
    for (size_t i = 0; i < natoms; i++)
    {
        for (size_t j = 0; j < natoms; j++)
        {
            if (j != i)
            {
                dist = get_dist(i, j);

                if (dist <= cutoff)
                {
                    coordNums[i]++;
                }
            }
        }
    }
}



// Get vector between atom i and j
xyz Trajectory::get_vec(int i, int j) // use get_dist to calculate MIC
{
    xyz ij;                // vector ij
    
    ij.x = coords[j].x - coords[i].x;
    ij.y = coords[j].y - coords[i].y;
    ij.z = coords[j].z - coords[i].z;

    // Minimum image convention    
    ij.x -= boxdims.x * round(ij.x / boxdims.x);
    ij.y -= boxdims.y * round(ij.y / boxdims.y);
    ij.z -= boxdims.z * round(ij.z / boxdims.z);

    return ij;
}


// Get vector between point i and j // for debugging
xyz Trajectory::get_vec(const xyz& i, const xyz& j)
{
    xyz ij;           // vector ij
    ij.x = j.x - i.x;
    ij.y = j.y - i.y;
    ij.z = j.z - i.z; 

    return ij;
}


// Calculate angle ijk between vector ij and ik (unit: degree)
double Trajectory::get_angle(const xyz& ij, const xyz& ik)
{
    double angle;     // angle between vector ij and ik 
    double rij, rik;  // L2 norm (modulus) of rij and rik

    rij = l2norm(ij);
    rik = l2norm(ik);
    
    // Calculate angle with arc cosine function (unit: rad)    
    angle = acos(dot(ij, ik) / (rij * rik)); // This value is always non-negative
    angle = angle * 180.0 / M_PI; // rad to degree
    
    return angle;
}



// Calculate the azimuth between vector il and plane determined by vector ij and ik (unit: degree) 
double Trajectory::get_azimuth(const xyz& ij, const xyz& ik, const xyz& il)
{   
    // Normal vector of plane determined by ij and ik
    // NOTE THAT IJ AND IK SHOULD NEVER BE THE SAME OR THE OPPOSITE VECTOR!!!!!
    xyz normVec = cross(ij, ik); 

    // Azimuth between vector il and plane deteremined by vector ij and ik 
    double azimuth;



    // Calculate the angle between il and the projection of il in plane determined by ij and ik


    // Check dot product of il and normVec
    if (dot(il, vec_add(ij, ik)) >= 0) // when azimuth <= 90
    {
        if (dot(il, normVec) >= 0) // angle between norm vec and il <= 90
            azimuth = 90.0 - get_angle(il, normVec); 
        else                       // angle between norm vec and il > 90 
            azimuth = get_angle(il, normVec) - 90.0;

    }
    else  // when the azimuth is > 90
    {
        if (dot(il, normVec) >= 0) // angle between norm vec and il <= 90
            azimuth = 90.0 + get_angle(il, normVec); 
        else                       // angle between norm vec and il > 90  
            azimuth = 270.0 - get_angle(il, normVec);
    }

    return azimuth;    




    // cout << fixed << setprecision(4) << "Normal vector of plane 0102: " << normVec.x << " " << normVec.y << " " << normVec.z << endl;  // VALUE PRINTED LOOKS RIGHT!

    // cout << fixed << setprecision(4) << "Vector of 03: " << il.x << " " << il.y << " " << il.z << endl; 
        
    // cout << fixed << setprecision(4) << "Angle between 03 and normal vector of plane 0102: " << get_angle(il, normVec) << endl; // VALUE PRINTED IS NOT RIGHT!!





    // // Projection of il to plane determined by ij and ik
    // xyz proj_il_ijik;

    // // Calculate the projection of il in plane determined by ij and ik
    // double norm_proj_il_ijik = 1.0 / (l2norm(normVec) * l2norm(normVec));  //norm of projection

    // xyz unitVec_proj_il_ijik = cross(normVec, cross(il, normVec)); // projection unit vector

    // proj_il_ijik = scal_prod(unitVec_proj_il_ijik, norm_proj_il_ijik);

    // cout << "Projection vector: " << proj_il_ijik.x << " " << proj_il_ijik.y << " " << proj_il_ijik.z << endl;

    // azimuth = get_angle(il, proj_il_ijik);
}



// Get neighbor list of each atom in the frame
void Trajectory::get_neighlist(double cutoff)
{
    neighlist.clear();
    neighlist.resize(natoms); 

    double dist;
    
    for (size_t i = 0; i < natoms; i++)
    {
        for (size_t j = 0; j < natoms; j++)
        {
            if (j != i)
            {
                dist = get_dist(i, j);

                if (dist <= cutoff)
                {
                    neighlist[i].push_back(j);
                }
            }
        }
    }
}



// Get sp3 order parameter of each atom in the frame 
void Trajectory::get_sp3OP(double del_theta) 
// del_theta is the tolerance angle 
{
    sp3OP.clear();
    sp3OP.resize(natoms, 0.0);

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in sp3OP
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of sp3OP

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp3OP
    double exp_k = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)> 


    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {
       
        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  
        
        //for debugging
        //cout << endl;
        //cout << "Neighbor list of atom " << i << ": ";
        //for (size_t ii = 0; ii < nneighbor; ii++)
        //{
        //    cout << neighlist[i][ii] << " ";
        //}
        //cout << endl;


        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp3OP[i] = 0.0;
          continue;
        }

        // // for debugging
        // cout << endl;
        // cout << "############################" << endl;
        // cout << "Atom number: " << i << endl;
        // cout << "Number of neighbors of atom " << i << ": " << nneighbor << endl;

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);
        
        // for debugging
        // cout << fixed << setprecision(4) << "Prefactor of atom i: " << i << ": " << prefactor << endl;


        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                cossq_exp_l = 0.0;

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for debugging 
                //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][k] << ": " << angle_ij_ik << endl;



                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    // for debugging
                    //cout << "i j k l: " << i << " " << neighlist[i][j] << " " << neighlist[i][k] << " " << neighlist[i][l] << endl;




                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    // for debugging
                    //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][l] << ": " << angle_ij_il << endl;




                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    
                    //for debugging
                    //cout << fixed << setprecision(4) << "Azimuth_" << i << neighlist[i][l] << "_" << i << neighlist[i][j] << i << neighlist[i][k] << ": " << azimuth_il_ijik << endl;



                    // populate cossq_exp_l (inner summation)
                    cossq_exp_l += pow(cos(M_PI / 180.0 * 1.5 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(del_theta, 2))); 

                    // // for debugging
                    // cout << "cossq: " << pow(cos(1.5 * azimuth_il_ijik), 2) << endl;
                    // cout << "exp_l: " << exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

                }


                // populate exp_k (outer summation)
                exp_k += exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(del_theta, 2))) * cossq_exp_l;

                // // for debugging
                // cout << "exp_k: " << exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

            }

            // // calculate and populate (for averaging) sp3OP of atom i
            // cout << fixed << setprecision(4) << "sp3OP: " << prefactor * exp_k << endl;
            // cout << endl;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp3OP[i] += exp_k;
        }

        // calculate sp3OP of atom i
        sp3OP[i] *= prefactor;

        // for debugging
        //cout << "sp3OP of atom " << i << ": " << sp3OP[i] << endl;
    }    

}



// Get sp2 order parameter of each atom in the frame 
void Trajectory::get_sp2OP(double del_theta) 
// del_theta is the tolerance angle 
{
    sp2OP.clear();
    sp2OP.resize(natoms, 0.0);

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in sp2OP
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of sp3OP

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp2OP
    double exp_k = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)> 


    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {

        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  

        //for debugging
        //cout << endl;
        //cout << "Neighbor list of atom " << i << ": ";
        //for (size_t ii = 0; ii < nneighbor; ii++)
        //{
        //    cout << neighlist[i][ii] << " ";
        //}
        //cout << endl;


        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp2OP[i] = 0.0;
          continue;
        }   

        // // for debugging
        // cout << endl;
        // cout << "############################" << endl;
        // cout << "Atom number: " << i << endl;
        // cout << "Number of neighbors of atom " << i << ": " << nneighbor << endl;

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);
        
        // for debugging
        // cout << fixed << setprecision(4) << "Prefactor of atom i: " << i << ": " << prefactor << endl;


        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                cossq_exp_l = 0.0;

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for debugging 
                //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][k] << ": " << angle_ij_ik << endl;



                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    // for debugging
                    //cout << "i j k l: " << i << " " << neighlist[i][j] << " " << neighlist[i][k] << " " << neighlist[i][l] << endl;




                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    // for debugging
                    //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][l] << ": " << angle_ij_il << endl;




                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    
                    //for debugging
                    //cout << fixed << setprecision(4) << "Azimuth_" << i << "-" << neighlist[i][l] << "_" << i << "-" << neighlist[i][j] << "-" << i << "-" << neighlist[i][k] << ": " << azimuth_il_ijik << endl;



                    // populate cossq_exp_l (inner summation)
                    cossq_exp_l += pow(cos(M_PI / 180.0 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 120.0), 2) / (2 * pow(del_theta, 2))); 

                    // // populate cossq_exp_l (inner summation)
                    // cossq_exp_l += pow(cos(azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 120.0), 2) / (2 * pow(del_theta, 2))); 



                    // // for debugging
                    // cout << "cossq: " << pow(cos(1.5 * azimuth_il_ijik), 2) << endl;
                    // cout << "exp_l: " << exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

                }


                // populate exp_k (outer summation)
                exp_k += exp(- pow((angle_ij_ik - 120.0), 2) / (2 * pow(del_theta, 2))) * cossq_exp_l;

                // // for debugging
                // cout << "exp_k: " << exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

            }

            // // calculate and populate (for averaging) sp3OP of atom i
            // cout << fixed << setprecision(4) << "sp3OP: " << prefactor * exp_k << endl;
            // cout << endl;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp2OP[i] += exp_k;
        }

        // calculate sp3OP of atom i
        sp2OP[i] *= prefactor;

        // for debugging
        //cout << "sp2OP of atom " << i << ": " << sp3OP[i] << endl;
    }    

}


// Get sp2 and sp3 order parameter of each atom in the frame 
void Trajectory::get_OP(double sp3_del_theta, double sp2_del_theta) 
// del_theta is the tolerance angle 
{
    sp2OP.clear();
    sp3OP.clear();
    sp2OP.resize(natoms, 0.0);
    sp3OP.resize(natoms, 0.0);

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in sp2OP
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of sp3OP

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp2OP
    double exp_k_2 = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l_2 = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)>

    // Terms in sp3OP
    double exp_k_3 = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l_3 = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)>  


    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {

        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  

        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp2OP[i] = 0.0;
          sp3OP[i] = 0.0;
          continue;
        }   

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);

        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k_2 = 0.0;
            exp_k_3 = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                cossq_exp_l_2 = 0.0;
                cossq_exp_l_3 = 0.0; 

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    
                    // populate cossq_exp_l_2 (inner summation)
                    cossq_exp_l_2 += pow(cos(M_PI / 180.0 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 120.0), 2) / (2 * pow(sp2_del_theta, 2)));

                    // populate cossq_exp_l_3 (inner summation)
                    cossq_exp_l_3 += pow(cos(M_PI / 180.0 * 1.5 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(sp3_del_theta, 2)));  

                }

                // populate exp_k_2 (outer summation)
                exp_k_2 += exp(- pow((angle_ij_ik - 120.0), 2) / (2 * pow(sp2_del_theta, 2))) * cossq_exp_l_2;

                // populate exp_k_3 (outer summation)
                exp_k_3 += exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(sp3_del_theta, 2))) * cossq_exp_l_3;

            }

            // sum up the terms of sp2OP of atom i in each loop of j
            sp2OP[i] += exp_k_2;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp3OP[i] += exp_k_3;
        }

        // calculate sp2OP of atom i
        sp2OP[i] *= prefactor;

        // calculate sp3OP of atom i
        sp3OP[i] *= prefactor;
    }
}