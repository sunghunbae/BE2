#include <vector>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "lib.h"

using namespace std;

void init_fyc_lib (vector < struct lib > & fyc) 
{
    struct lib r;
    r.f.clear();
    r.y.clear(); 
    r.fP.clear();
    r.yP.clear(); 
    r.c1.clear();
    r.c2.clear();
    r.c3.clear();
    r.c4.clear();

    r.seq = "GLY"; fyc.push_back (r); // 0
    r.seq = "ALA"; fyc.push_back (r); // 1
    r.seq = "VAL"; fyc.push_back (r); // 2
    r.seq = "LEU"; fyc.push_back (r); // 3
    r.seq = "ILE"; fyc.push_back (r); // 4
    r.seq = "MET"; fyc.push_back (r); // 5
    r.seq = "TYR"; fyc.push_back (r); // 6
    r.seq = "ASP"; fyc.push_back (r); // 7
    r.seq = "ASN"; fyc.push_back (r); // 8
    r.seq = "GLU"; fyc.push_back (r); // 9
    r.seq = "GLN"; fyc.push_back (r); // 10
    r.seq = "CYS"; fyc.push_back (r); // 11
    r.seq = "SER"; fyc.push_back (r); // 12
    r.seq = "THR"; fyc.push_back (r); // 13
    r.seq = "TRP"; fyc.push_back (r); // 14
    r.seq = "PHE"; fyc.push_back (r); // 15
    r.seq = "ARG"; fyc.push_back (r); // 16
    r.seq = "HIS"; fyc.push_back (r); // 17
    r.seq = "LYS"; fyc.push_back (r); // 18
    r.seq = "PRO"; fyc.push_back (r); // 19
}

size_t read_FYCQR (const char *file, 
    vector < struct lib > & fyc, vector < struct ff > & qr) 
{
    vector < string > c;
    struct ff amber99;
    string buffer,fyPRO;
    int i,n=0;

    ifstream lib (file);
    if (lib.is_open()) {
	while (!lib.eof()) {
	    getline (lib,buffer);
	    if (parse(buffer.c_str(), "\t ", c) != 0) {

		// phi-psi library
		if (c[0] == "FY") {
		   for (i=0;i<fyc.size();i++) {
			fyPRO = fyc[i].seq + "PRO";
			if (c[1] == fyc[i].seq) {
			    n++;
			    fyc[i].f.push_back(c[2]);
			    fyc[i].y.push_back(c[3]);
			    break;
			    } 
			if (c[1] == fyPRO) {
			    n++;
			    fyc[i].fP.push_back(c[2]);
			    fyc[i].yP.push_back(c[3]);
			    break;
			    }
			} // for
		    } // phi-psi library

		// chi library
		if (c[0] == "CX") {
		    for (i=0;i<fyc.size();i++) {
			if (c[1] == fyc[i].seq) {
			    n++;
			    fyc[i].c1.push_back(c[2]);
			    fyc[i].c2.push_back(c[3]);
			    fyc[i].c3.push_back(c[4]);
			    fyc[i].c4.push_back(c[5]);
			    break;
			    }
			}
		    } // chi library

		// partial charge and VDW radius (amber99.lib)
		if (c[0] == "QR") {
		    amber99.residue = c[1];
		    amber99.atom = c[2];
		    sscanf(c[3].c_str(),"%lf",&amber99.charge);
		    sscanf(c[4].c_str(),"%lf",&amber99.radius);
		    qr.push_back (amber99);
		    } // partial charge and VDW radius

		} // if not comment
	    } // while

	lib.close();

	} // if opened
    else {
	cout << "error: cannot open " << file << endl; 
	exit(1);
	}

    return n;
}

size_t parse (const char *line, const char *delimit, vector <string> &c)
{
  char work[strlen(line)],*token;
  bool comment = false;
  int i;

  strcpy (work, line);
  c.resize(0);

  token = strtok (work, delimit);
  if (token == NULL) return 0;

  for (i=0;i<strlen(token);i++)
    if (token[i] == '#') {
      comment = true;
      break;

      }

  if (!comment)
    c.push_back (token);
  else
    return 0;

  while (token != NULL) {
    token = strtok (NULL, delimit);
    if (token != NULL) {
      for (i=0;i<strlen(token);i++)
        if (token[i] == '#') {
          comment = true;
          break;
          }
      if (comment)
        break;
      else
        c.push_back (token);
      }
    }
  return c.size();
}
