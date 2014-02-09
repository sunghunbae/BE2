#define TITLE "Ensemble Generator written by Sung-Hun Bae"
#define SIGNATURE "June-30-2009"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "lib.h"
#include "pdb.h"
#include "clash.h"
#include "rotate.h"

#define MAXLEN 1024

void usage ();
