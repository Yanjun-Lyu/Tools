#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "helpers.hpp"

using namespace std;

bool get_next_line(ifstream & instream, string & line)
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
// The program will not be terminated even if encountering an error reading the line.
// Use it when you want to continue execution when encountering an error (for uncritical input reading)
{
    // getline(instream, line);
    // if (! instream.good())
    //     return false;
    // else
    //     return true;

    if (getline(instream, line)) // Attempt to read a line
        return true; // If successful, return true
    else if (!instream.eof()) // If getline() failed but EOF flag is not set
    {
        cout << "Encountered unexpected error while reading file." << endl;
        exit(1);
    }
    return false; // If getline() failed and EOF flag is set, return false
}


string get_next_line(ifstream & instream)
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
// The program will be terminated if encountering an error reading the line.
// Use it when you want to terminate the program when encountering an error. (For critical input reading)
{
    // string line;
    
    // getline(instream, line) ;
    // if (! instream.good())
    // {
    //     cout << "Encountered unexpected error while reading file." << endl;
    //     exit(1);
    //     return line;
    // }
    // else
    //     return line;

    string line;
    
    if (getline(instream, line)) // Attempt to read a line
        return line; // If successful, return the line
    else if (!instream.eof()) // If getline() failed but EOF flag is not set
    {
        cout << "Encountered unexpected error while reading file." << endl;
        exit(1);
    }
    return ""; // If getline() failed and EOF flag is set, return an empty string

}


int tokenize(string line, vector<string>& tokens)
// Break a line up into tokens based on space separators. Returns the number of tokens parsed.
{
    string buf;
  
    stringstream stream_parser;

    int pos = line.find('\n');   // Strip off terminal new line.
    if ( pos != string::npos ) 
        line.erase(pos, 1);

    stream_parser.str(line);
	 
    tokens.clear();

    while ( stream_parser >> buf ) 
        tokens.push_back(buf);

    return(tokens.size() );
}


string get_token(string line, int field)
// Returns a specific token from a line
{
    vector<string> tokens;
    tokenize(line, tokens);
    
    return tokens[field];
}



//////////////////////////////////////////////////////////////////////
// 3D Vector calculation
//////////////////////////////////////////////////////////////////////
// vecs are passed by reference for efficiency
// const keyword is used to make sure they are not modified in the function
// same implementation scheme is adapted in all vector calculation functions

// Vector add
xyz vec_add(const xyz& vec1, const xyz& vec2)
{
    xyz vec;
    vec.x = vec1.x + vec2.x; 
    vec.y = vec1.y + vec2.y;
    vec.z = vec1.z + vec2.z;
    return vec; 
}


// Vector minus (vec1 - vec2)
xyz vec_minus(const xyz& vec1, const xyz& vec2)
{
    xyz vec;
    vec.x = vec1.x - vec2.x; 
    vec.y = vec1.y - vec2.y;
    vec.z = vec1.z - vec2.z;
    return vec; 
}


// Scalar product of a vector and a number
xyz scal_prod(const xyz& vec, double n)
{
    xyz prod;
    prod.x *= n; 
    prod.y *= n;
    prod.z *= n;
    return prod; 
}


// L2 norm (vector mod)
double l2norm(const xyz& vec)
{
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}


// Dot product
double dot(const xyz& vec1, const xyz& vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}


// Cross product
xyz cross(const xyz& vec1, const xyz& vec2) 
{
    xyz result;
    result.x = vec1.y * vec2.z - vec1.z * vec2.y;
    result.y = vec1.z * vec2.x - vec1.x * vec2.z;
    result.z = vec1.x * vec2.y - vec1.y * vec2.x;
    if (result.x == 0.0 && result.y == 0.0 && result.z == 0.0)
    {
        cerr << "Same/Opposite vectors are used to calculate cross product, which causes nan when calculate azimuths! " << endl;
    } 
    return result;
}

