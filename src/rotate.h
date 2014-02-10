#ifndef __ROTATE__
#define __ROTATE__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "lib.h"
#include "pdb.h"
#include "topo.h"

typedef	struct {
	double x,y,z;
	} Vector;

int	read_rot_list (const char *file, vector < struct rot > &rlist);
void	write_rot_list (const char *prefix, const vector < struct rot > &rlist,
	const int idx, const int num_clash, const double Evdw);
void	rotate (vector <ATOM> &coor, vector <struct rot> &R, 
	vector <struct lib> &fyc);
void	dihedrot (vector <ATOM> &, const struct rot &);
void	get_sequence (const vector <ATOM> & coor,
        const char chainId, const unsigned int resSeq, string &);
void	get_vector (const vector <ATOM> & coor,
        char chainId, unsigned int resSeq, string p, Vector &a);

int	RotateAtom_ (ATOM &, Vector *v1, Vector *v2, double angle);
double	Dihedral_ (Vector *a, Vector *b, Vector *c, Vector *d);
void	VectorProduct_ (Vector *v3,Vector *v1,Vector *v2);
double	ScalarProduct_ (Vector *v1, Vector *v2); 
double	AbsoluteValue_ (Vector *v); 
int	ParallelPart_ (Vector *par_v, Vector *ref_v, Vector *input_v);

#endif
