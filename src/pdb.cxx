#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

using namespace std;

#include "pdb.h"
#include "lib.h"

size_t read_PDB (const char *filename, vector <ATOM> & coor)
{
  char line [128];
  string s,temp,temp_;
  size_t c;

  ifstream pdb_file (filename);
  if (pdb_file.fail()) {
    printf("error: cannot open %s\n",filename);
    exit (-1);
    }

  coor.resize(0);
  ATOM a;

  while (! pdb_file.eof()) {
    pdb_file.getline(line, sizeof(line));
    if (strlen(line) > 0) {
      s = string(line);
      if(s.find("ATOM") != string::npos || s.find("HETATM") != string::npos) {

        // PDB format version 2.3
	temp	    = s.substr(6,5);  
	sscanf(temp.c_str(),"%d",&a.serial);
        temp	    = s.substr(12,4);
	a.name	    = temp;

	while ((c=temp.find(" ",0)) != string::npos) 
	  temp.replace (c,1,""); // trim trailing spaces
	a.name_	    = temp;
	temp_	    = temp.substr(0,1);
	if (temp_ == "1" || temp_ == "2" || temp_ == "3")
	  a.name_   = temp.substr(1) + temp_ ;

	temp	    = s.substr(16,1); 
	sscanf(temp.c_str(),"%c",&a.altLoc);

        a.resName   = s.substr(17,3);
	a.resName_  = a.resName;
	// exceptions - amber99 names
	if (a.name_ == "HT1" || a.name_ == "HT2" || a.name_ == "HT3") {
	  a.resName_ = "N" + a.resName_;
	  if (a.name_ == "HT1") a.name_ = "H1";
	  if (a.name_ == "HT2") a.name_ = "H2";
	  if (a.name_ == "HT3") a.name_ = "H3";
	  }
	if (a.name_ == "OXT") {
	  a.resName_ = "C" + a.resName_;
	  }
	if (a.resName_ == "GLY" && a.name_ == "HA")
	  a.name_   = "HA2";
	if (a.resName_ == "HIS")
	  a.resName_ = "HIE";

        temp   	    = s.substr(21,1); 
	sscanf(temp.c_str(),"%c",&a.chainId);
        temp  	    = s.substr(22,4);
	sscanf(temp.c_str(),"%d",&a.resSeq);
	temp 	    = s.substr(26,1);
	sscanf(temp.c_str(),"%c",&a.iCode);
        temp	    = s.substr(30,24); 
	sscanf(temp.c_str(),"%lf %lf %lf", &a.x, &a.y, &a.z);
	temp        = s.substr(54,6);
	sscanf(temp.c_str(),"%lf", &a.occupancy);
	temp        = s.substr(60,6);
	sscanf(temp.c_str(),"%lf", &a.tempFactor);
	a.segID	    = s.substr(72,4);

	temp	    = s.substr(76,2);
	a.element   = temp;
	while ((c=temp.find(" ",0)) != string::npos) 
	  temp.replace (c,1,""); // trim trailing spaces
	a.element_  = temp;

	a.charge    = s.substr(78,2);
	
        coor.push_back(a);

        } // ATOM or HETATM
      } // if line is not empty
    } // until end of file

  return coor.size();
}

void write_PDB (const char *prefix, const vector <ATOM> &coor, 
    const int idx, const int num_clash, const double Evdw)
{
    int i;
    char filename [256];
    sprintf(filename, "%s_%04d.pdb", prefix, idx);
    FILE * fileO = fopen(filename,"w");

    // REMARK
    if (num_clash >= 0)
	fprintf (fileO, "REMARK clash %d Evdw %.2f\n",num_clash,Evdw);

    // PDB version 2.3
    for (i=0; i<coor.size(); i++) {
	fprintf(fileO, "ATOM  %5d %-4s%c%3s ",
	    coor[i].serial, coor[i].name.c_str(), coor[i].altLoc, 
	    coor[i].resName.c_str());
	fprintf(fileO, "%c%4d%c   ",
	    coor[i].chainId, coor[i].resSeq, coor[i].iCode);
	fprintf(fileO, "%8.3f%8.3f%8.3f",
	    coor[i].x, coor[i].y, coor[i].z);
	fprintf(fileO, "%6.2f%6.2f      ",
	    coor[i].occupancy, coor[i].tempFactor);
	fprintf(fileO, "%4s%2s%2s\n",
	    coor[i].segID.c_str(), coor[i].element.c_str(), 
	    coor[i].charge.c_str());
	}
    fclose(fileO);
}

void write_PQR (const char *prefix, const vector <ATOM> &coor, const int idx)
{
  int i;
  char filename [256];

  sprintf(filename, "%s_%04d.pdb", prefix, idx);

  FILE * fileO = fopen(filename,"w");

  /* output (according to PDB version 2.3) */
  for (i=0; i<coor.size(); i++) {
    fprintf(fileO, "ATOM  %5d %-4s%c%3s ",
      coor[i].serial, coor[i].name.c_str(), coor[i].altLoc, 
      coor[i].resName.c_str());
    fprintf(fileO, "%c%4d%c   ",
      coor[i].chainId, coor[i].resSeq, coor[i].iCode);
    fprintf(fileO, "%8.3f%8.3f%8.3f",
      coor[i].x, coor[i].y, coor[i].z);
    fprintf(fileO, "%6.2f%6.2f      ",
      coor[i].charge_, coor[i].radius_);
    fprintf(fileO, "%4s%2s%2s\n",
      coor[i].segID.c_str(), coor[i].element.c_str(), 
      coor[i].charge.c_str());
    }

  fclose(fileO);
}
