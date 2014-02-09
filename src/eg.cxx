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

#define SIGNATURE "Sung-Hun Bae 2010.07"
#define MAXLEN 1024

void usage ();

using namespace std;

unsigned long int seed;
gsl_rng * rng;

double vdw_overlap = 1.6;
double hbond_limit = 1.6;

int main (int argc, char *argv[])
{
    // setup random number generator
    extern unsigned long int seed;
    extern gsl_rng * rng;
    const gsl_rng_type * T = gsl_rng_ranlux389;
    rng = gsl_rng_alloc (T);

    // default seed
    seed = time(0);

    char file_lib [MAXLEN];
    char file_pdb [MAXLEN];
    char file_rot [MAXLEN];
    char o_prefix [MAXLEN];

    int i,n,e;
    int o_num_begin,o_num_end;
    int max_allowed_clash, clash_0, clash_1;
    double max_allowed_Evdw, Evdw;

    vector < struct lib > fyc;
    vector < struct ff > qr;
    vector < struct rot > rlist;
    vector < struct exl > xlist;
    vector < ATOM > coor;
    vector < ATOM > coor_ini;

    init_fyc_lib (fyc);

    /* default setting */
    strcpy (file_lib, "./fycqr.lib");
    strcpy (file_rot, "./rlist");
    strcpy (o_prefix, "e");
    max_allowed_clash	= 0;
    max_allowed_Evdw	= 10000.;
    o_num_begin		= 1;
    o_num_end		= 1;
    bool _quiet		= false;
    bool _import_lib	= true;
    bool _test_clash_0	= true;
    bool _test_clash_1	= true;
    bool _save_rlist	= false;
    bool _save_pdb	= false;
    bool _rebuild	= false;

    if (argc < 2) usage();

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-h") == 0) usage();
	if (strcmp(argv[i],"-i") == 0) strcpy(file_pdb, argv[++i]);
	if (strcmp(argv[i],"-r") == 0) strcpy(file_rot, argv[++i]);
	if (strcmp(argv[i],"-l") == 0) strcpy(file_lib, argv[++i]);
	if (strcmp(argv[i],"-d") == 0) {
	    // make a list for exclusin of the clash test
	    struct exl x;
	    sscanf(argv[++i], "%c", &x.chain);
	    sscanf(argv[++i], "%d", &x.nter);
	    sscanf(argv[++i], "%d", &x.cter);
	    xlist.push_back (x); 
	    }
	if (strcmp(argv[i],"-maxc") == 0) {
	    sscanf(argv[++i], "%d", &max_allowed_clash);
	    }
	if (strcmp(argv[i],"-maxE") == 0) {
	    sscanf(argv[++i], "%lf", &max_allowed_Evdw);
	    }
	if (strcmp(argv[i],"-o") == 0) {
	    strcpy(o_prefix, argv[++i]);
	    sscanf(argv[++i], "%d", &o_num_begin);
	    sscanf(argv[++i], "%d", &o_num_end);
	    }
	if (strcmp(argv[i],"-rebuild") == 0) {
	    _rebuild	    = true;
	    _import_lib	    = false;
	    _test_clash_0   = false;
	    _test_clash_1   = false;
	    _save_pdb	    = true;
	    _save_rlist	    = false;
	    clash_1	    = -1;
	    Evdw	    = -1.;
	    }
	if (strcmp(argv[i],"-seed") == 0) sscanf(argv[++i], "%lu", &seed);
	if (strcmp(argv[i],"-pdb") == 0) _save_pdb = true;
	if (strcmp(argv[i],"-rlist") == 0) _save_rlist = true;
	if (strcmp(argv[i],"-quiet") == 0) _quiet = true;
	i++;
	}

    // setup random number generator
    gsl_rng_set (rng, seed);

    // import phi/psi angles and van der Waals radius from library
    if (_import_lib) {
	if (!_quiet) {
	    clog << "# library: ";
	    clog << string (file_lib);
	    clog << " (" << read_FYCQR (file_lib, fyc, qr) << ")";
	    clog << endl;
	    }
	else
	    read_FYCQR (file_lib, fyc, qr);
	}

    // read PDB file
    if (!_quiet) {
	clog << "# input PDB: " << string (file_pdb);
	clog << " (" << read_PDB (file_pdb, coor) << ")";
	clog << endl;
	}
    else
	read_PDB (file_pdb, coor);

    // read rotation list file
    if (!_quiet) {
	clog << "# rotation list: " << string (file_rot);
	clog << " (" << read_rot_list (file_rot, rlist) << ")";
	clog << endl;
	}
    else
	read_rot_list (file_rot, rlist);

    // rebuild pdb from rotation list output (i.e. prefix_####.rot)
    if (_rebuild) {
	string str = string(file_rot);
	size_t found1 = str.find_last_of ("_");
	size_t found2 = str.find_last_of (".");
	string prefix = str.substr(0,found1);
	string number = str.substr(found1+1,found2);
	sscanf (prefix.c_str(),"%s",o_prefix);
	sscanf (number.c_str(),"%d",&o_num_begin);
	o_num_end = o_num_begin;
	}

    // set partial charge (Q) and van der Waals radius (R)
    if (_test_clash_0 || _test_clash_1) 
	set_QR (coor, qr);

    // keep the original coordinates
    coor_ini = coor; 

    // check initial steric clash
    if (_test_clash_0) {
	if (!_quiet) {
	    clog << "# testing steric clash of the input PDB coordinates"; 
	    clog << endl;
	    clash_0 = test_steric_clash (coor_ini, xlist, Evdw,true);
	    }
	else
	    clash_0 = test_steric_clash (coor_ini, xlist, Evdw,false);
	}

    // loop for ensemble generation
    for(e=o_num_begin;e<=o_num_end;e++) {

	// reset coordinates
	coor = coor_ini;

	// apply rotation list & check steric clash
	do {
	    rotate (coor, rlist, fyc);
	    if (_test_clash_1) 
		clash_1 = test_steric_clash (coor, xlist, Evdw, false);
	    } while (_test_clash_1 && 
		(clash_1 > max_allowed_clash || Evdw > max_allowed_Evdw));

	// clash test results
	if (!_quiet && _test_clash_1) {
	    clog << e;
	    clog << " clash= " << clash_1;
	    clog << fixed << setprecision(2);
	    clog << " Evdw= " << Evdw;
	    clog << "  ";
	    }

	// output rotation list (.rot)
	if (_save_rlist) {
	    write_rot_list (o_prefix, rlist, e, clash_1, Evdw);
	    if (!_quiet) {
		clog << o_prefix <<"_"<<setw(4)<<setfill('0') << e;
		clog << ".rot ";
		}
	    }

	// output pdb (.pdb)
	if (_save_pdb) {
	    write_PDB (o_prefix, coor, e, clash_1, Evdw);
	    if (!_quiet) {
		clog << o_prefix <<"_"<<setw(4)<<setfill('0') << e;
		clog << ".pdb ";
		}
	    }

	if (!_quiet && (_test_clash_1 || _save_pdb || _save_rlist))
	    clog << endl;

	} // e, ensemble

    gsl_rng_free (rng);
}

void usage ()
{
  printf("\t\t\t\t\t\t\t" SIGNATURE "\n\n");
  printf("  Ensemble Generator (coil library/dihedral angle rotations)\n");
  printf("\n");
  printf("  Usage: eg [options]\n");
  printf("  options:\n");
  printf("    -i #        read <i>nput PDB file\n");
  printf("    -d # # #    rigid <d>omain [chainId begin end]\n");
  printf("    -r #        read <r>otation list (default= rlist)\n");
  printf("    -l #        import <l>ibrary (default= fycqr.lib)\n");
  printf("    -o # # #    <o>utput [prefix begin end]\n");
  printf("    -maxc #     <max>imum allowed <c>lashes (default= 0)\n");
  printf("    -maxE #     <max>imum allowed VDW <e>nergy (default= 10000)\n");
  printf("    -seed #     <seed> for random number generation\n");
  printf("    -pdb        save output as <PDB> (.pdb)\n");
  printf("    -rlist      save output as <rlist> (.rot)\n");
  printf("    -rebuild    re<b>uild from rlist\n");
  printf("    -quiet      be <quiet> do not print anything\n");
  printf("\n");
  printf("  Example: eg -i a.pdb -d A 8 62 -d A 90 145 -o x 1 10\n");
  printf("\n");
  exit(1);
}
