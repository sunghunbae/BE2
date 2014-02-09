#ifndef __CLASH__
#define __CLASH__

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "pdb.h"
#include "lib.h"

void	set_QR (vector <ATOM> &coor, vector <struct ff> &amber99);
int	test_steric_clash (vector <ATOM> &coor, vector <struct exl> &, 
	double &, bool);
bool	is_excluded (ATOM &i, ATOM &j, vector <struct exl> &xlist); 
bool	cov_bonds_3 (ATOM &i, ATOM &j);

#endif
