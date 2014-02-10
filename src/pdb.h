#ifndef __PDB__
#define __PDB__

#include <string>
#include <vector>

using namespace std;

//
// PDB format version 2.3
//

typedef struct {
    unsigned int  serial, resSeq;
    string  name, resName, segID, element, charge;
    char    altLoc, chainId, iCode;
    double  x, y, z, occupancy, tempFactor;

    string  name_;  // without trailing spaces
    string  element_; // without trailing spaces
    string  resName_;
    string  atomtype_;
    double  radius_;
    double  charge_;
    } ATOM;

unsigned int	read_PDB (const char *file, vector <ATOM> & coor);
void	write_PDB (const char *prefix, const vector <ATOM> &coor, 
	    const int idx, const int num_clash, const double Evdw);
void	write_PQR (const char *prefix, const vector <ATOM> &coor, 
	    const int idx);

#endif
