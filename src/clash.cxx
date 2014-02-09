#include "clash.h"

using namespace std;

// 
//  clash_limit = sum of van der Waals radii - overlap
//  overalp     = 0.6-1.6 A (default=1.2)
//  clash, if distance is less than clash_limit
//
//  ignore atom pairs of fewer than 3 bonds apart (4 bonds apart for H)
//

void set_QR (vector <ATOM> &coor, vector <struct ff> &amber99)
{
    int i,f;
    for (i=0; i<coor.size(); i++) {
        for (f=0; f<amber99.size(); f++) {
            if (coor[i].resName_ == amber99[f].residue &&
                coor[i].name_ == amber99[f].atom ) {
                coor[i].charge_ = amber99[f].charge;
                coor[i].radius_ = amber99[f].radius;
                break;
                } // if
            } // f
        } // i
}

int test_steric_clash (vector <ATOM> &coor, 
  vector <struct exl> &xlist, double &Evdw, bool verb)
{
    extern double vdw_overlap;
    extern double hbond_limit;
    double rij, vdwij, clash_limit;
    int i,j,clash_count = 0;
  
    Evdw = 0.0;

    for (i=0;i<coor.size();i++) {
	for (j=i+1;j<coor.size();j++) {
	    if (!((coor[i].chainId == coor[j].chainId) &&
		    (coor[i].resSeq == coor[j].resSeq)) &&
		(coor[i].x - coor[j].x) < 3. &&
		(coor[i].y - coor[j].y) < 3. &&
		(coor[i].z - coor[j].z) < 3. &&
		! is_excluded (coor[i],coor[j],xlist)) {

		vdwij = coor[i].radius_ + coor[j].radius_;

		/* exception: within 3 covalent bonds */
		if (cov_bonds_3 (coor[i],coor[j]))
		    clash_limit = 0.0;
		/* allow hydrogen bonding */
		else if ((coor[i].element_ == "H" && 
		    (coor[j].element_ == "N" || coor[j].element_ == "O")) ||
		    (coor[j].element_ == "H" && 
		    (coor[i].element_ == "N" || coor[i].element_ == "O")))
		    clash_limit = hbond_limit;
		else
		    clash_limit = vdwij - vdw_overlap;

		rij =  (coor[i].x - coor[j].x)*(coor[i].x - coor[j].x);
		rij += (coor[i].y - coor[j].y)*(coor[i].y - coor[j].y);
		rij += (coor[i].z - coor[j].z)*(coor[i].z - coor[j].z);
		rij = sqrt (rij);

		if (rij > vdwij) 
		    Evdw += 0.;
		else if (rij <= vdwij && rij >= 0.7*vdwij)
		    Evdw += -57.273*(1.0 - rij/(0.85*vdwij));
		else if (rij < 0.7*vdwij)
		    Evdw += 10.;

		if (rij < clash_limit) {
		    clash_count++;
		    if (verb) {
			clog << "#\t";
			clog << coor[i].resSeq << " ";
			clog << coor[i].name << " - ";
			clog << coor[j].resSeq << " ";
			clog << coor[j].name << " : ";
			clog << fixed << setprecision(2);
			clog << " rij " << rij; 
			clog << " limit " << clash_limit;
			clog << " vdwij " << vdwij;
			clog << endl;
			} // verb
		    } // clash

		} // clash test
	    } // j
	} // i

    if (verb) {
	clog << "#\tclash= " << clash_count;
	clog << fixed << setprecision(2);
	clog << " Evdw= " << Evdw << endl;
	}

    return clash_count;
}

//
// exclude rigid-rigid pairs from clash test
//
bool is_excluded (ATOM &i, ATOM &j, vector <struct exl> &xlist) 
{
  for (int k=0;k<xlist.size();k++) {
    if (i.chainId == xlist[k].chain && 
	i.resSeq  >= xlist[k].nter &&
	i.resSeq  <= xlist[k].cter &&
        j.chainId == xlist[k].chain && 
	j.resSeq  >= xlist[k].nter &&
	j.resSeq  <= xlist[k].cter)
      return true;
    }
  return false;
}


//
// return true if atoms i and j are within 3- or 4-(if hydrogren)
// covalent bonds distance
//
bool cov_bonds_3 (ATOM &i, ATOM &j)
{
    if (i.chainId != j.chainId) 
	return false;

    if ((i.resName == "CYS" || i.resName == "CYX") && i.name_ == "SG" &&
	(j.resName == "CYS" || j.resName == "CYX") && j.name_ == "SG" )
	return true;

    if (i.resSeq + 1 == j.resSeq) {
	if (i.name_ == "N"  && 
	    j.name_ == "N" || j.name_ == "HN" || j.name_ == "H" )
	    return true;
	if (i.name_ == "CA" && 
	    j.name_ == "N"  || j.name_ == "HN" || j.name_ == "H"  ||
	    j.name_ == "CA" || j.name_ == "HA" || 
	    (j.resName == "PRO" && j.name_ == "CD") )
	    return true;
	if (i.name_ == "CB" && 
	    j.name_ == "N" || j.name_ == "HN" || j.name_ == "H" )
	    return true;
	if (i.name_ == "C" && (
	    j.name_ == "N" || j.name_ == "HN" || j.name_ == "H" ||
	    j.name_ == "CA" || j.name_ == "HA" || 
	    j.name_ == "CB" || 
	    j.name_ == "HB1" || j.name_ == "HB2" || j.name_ == "HB" ||
	    j.name_ == "C"  || 
	    (j.resName == "PRO" && j.name_ == "CD") ||
	    (j.resName == "PRO" && j.name_ == "CG") ||
	    (j.resName == "PRO" && j.name_ == "2HG") ))
	    return true;
	if (i.name_ == "O" && (
	    j.name_ == "N"  || j.name_ == "HN" || j.name_ == "H" || 
	    j.name_ == "CA" || j.name_ == "HA" ||
	    (j.resName == "PRO" && j.name_ == "CD") ||
	    (j.resName == "PRO" && j.name_ == "CG") ||
	    (j.resName == "PRO" && j.name_ == "2HG") ))
	    return true;
	if (i.name_ == "CD" && 
	    j.resName == "PRO" && j.name_ == "N" )
	    return true;
	if (i.resName == "PRO" && (i.name_ == "C" || i.name_ == "O") &&
	    j.resName == "PRO" && j.name_ == "HD1" )
	    return true;
	}

    return false;
}
