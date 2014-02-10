#ifndef __LIB__
#define __LIB__

using namespace std;

struct lib { 
    string seq; 
    vector < string > f,y,fP,yP,c1,c2,c3,c4;
    };

struct rot {
    char chain;
    unsigned int resid;
    char rotation;
    bool library_search;
    string seq, dihed[6];
    };

struct exl {
  char chain;
  int nter;
  int cter;
  };

struct ff {
  string residue, atom;
  double charge, radius;
  };

void	init_fyc_lib (vector < struct lib > &);
unsigned int	read_FYCQR (const char *, vector < struct lib > &,
	vector < struct ff > &);
unsigned int	parse (const char *line, const char *delimit, vector <string> &c);

#endif
