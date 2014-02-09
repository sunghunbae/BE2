#define SIGNATURE "Sung-Hun Bae 2009.03"

#include <stdio.h>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

using namespace std;

typedef struct {
  double x,y,z;
  char c;
  int r;
  string name,resn;
  } Atom;

typedef struct {
  vector < vector < int > > atom;
  vector < char > chain;
  vector < int > resid;
  vector < vector < int > > anchor_element; // dR indices of anchor N/C ter
  vector < int > anchor_idx; // coor (pdb) index of anchor N/C ter
  vector < int > delta_r; // residue distance (|residue-residue|) to the anchor
  vector < double > d; // distance for velocity factor calc.
  } Attrib;

typedef struct {
  double ad1,ad2,ac,ab,bc,bd1,bd3,cd2,cd3,d1d2,a,sqa,sqb,sqc,sqd1,sqd2,sqd3;
  } scalar_product;

#define EPSABS     	0.0     /* integration, absolute error bound */
#define EPSREL     	1.0e-6  /* integration, relative error bound */
#define QAG_PREC_INTVL	50
#define QAG_RULE	GSL_INTEG_GAUSS21  /* options: 15,21,31,41,51,61 */
#define	MAXLENGTH	512
#define BOLTZMANN 	(1.3806503e-23)
#define C_TRANSLATION	(BOLTZMANN/(8*M_PI*1e-24))
#define C_ROTATION	(BOLTZMANN/(8*M_PI*1e-26))


#define _MSMS_FORMAT_	1
#define _GTS_FORMAT_	2
#define _OFF_FORMAT_	3

/* 
  Boltzmann constant kg m^2 K^-1 s^-2
  prefered unit for translational diffusion  = 10^-7 (cm)^2 s^-1
  prefered unit for rotational diffusion = 10^7 s^-1
*/

bool verb = false;	      /* default: short output */
bool gmap = false;	      /* map gboundary */
bool emap = false;	      /* map epsilon */
bool coff = false;	      /* GEOM colormap */

double temperature = 293.15;  /* 20 C, Kelvin */
double viscosity = 1.002;     /* water, cP. 1cP = 1e-3 kg m^-1 s^-1 */
double gboundary = 6.0;	      /* Amgstrom, 1A = 1e-10 m generous boundary */
double epsilon = 22.0;	      /* Amgstrom, 1A = 1e-10 m velocity decay coeff.*/
double map_g_min = 0.0;	      /* default map parameters */
double map_g_max = 30.0;      /* default map parameters */
double map_g_step = 2.0;      /* default map parameters */
double map_e_min = 1.0;	      /* default map parameters */
double map_e_max = 35.0;      /* default map parameters */
double map_e_step = 2.0;      /* default map parameters */

double assoc_threshold = 10.; /* default distance threshold for association */
int   algorithm = 1;	      /* default algorithm for velocity factor calc. */
char  Rchain_ = 'A';	      /* default chainId of rigid seg. */
int   RNter_ = 0;	      /* default N-terminus residue of rigid seg. */
int   RCter_ = 0;	      /* default C-terminus residue of rigid seg. */
vector < char > Rchain;
vector < int > RNter,RCter;


size_t read_msms_surface (const char *prefix, gsl_matrix *f);
size_t read_gts_surface (const char *prefix, gsl_matrix *f);
size_t read_geom_surface (const char *prefix, gsl_matrix *f);
size_t read_pdb_coor (const char *prefix, vector <Atom> & coor);
void write_geom_surface (const char *prefix, gsl_matrix *f,
    gsl_vector *vf, int color_scheme);
void calc_friction_tensor ( gsl_matrix *G_inv, gsl_matrix *C,
    gsl_vector *OR,gsl_vector *dS,gsl_vector *vf, 
    gsl_matrix *Ktt,gsl_matrix *Ktr, gsl_matrix *Krt,gsl_matrix *Krr);
void calc_diffusion_tensor (gsl_matrix *Ktt, gsl_matrix *Ktr,
    gsl_matrix *Krt, gsl_matrix *Krr, gsl_matrix *Dtt, 
    gsl_matrix *Dtr, gsl_matrix *Drt, gsl_matrix *Drr);
void display_matrix (
    gsl_matrix *Ktt, gsl_matrix *Ktr, gsl_matrix *Krt, gsl_matrix *Krr,
    gsl_matrix *Dtt, gsl_matrix *Dtr, gsl_matrix *Drt, gsl_matrix *Drr);
void calc_center_of_diffusion (gsl_matrix *Drr,gsl_matrix *Dtr,
    gsl_vector *OR);
void calc_surf_geometry (
    gsl_matrix *dSVert, gsl_matrix *C,
    gsl_matrix *dSSide, gsl_matrix *dSLength,gsl_vector *dS);
void calc_velocity_factor (gsl_matrix *C, Attrib &dSa, gsl_vector *vf);
void calc_attrib (gsl_matrix *C, gsl_matrix *RC, vector <Atom> &, Attrib &);
void calc_dimension (gsl_matrix *Vert); 
double dyt0 (double y, void * params);
double dyg0 (double y, void * params);
double dyyg0 (double y, void * params);
double dyyyg0 (double y, void * params);
double dyg1 (double y, void * params);
double dyyg1 (double y, void * params);
double dxxxg2 (double x, void * params);
void integ_err (size_t k,size_t j,scalar_product &); 

void block_inv   (gsl_matrix *G, gsl_matrix *G_inv);
void block_inv_2 (gsl_matrix *G, gsl_matrix *G_inv);
void block_inv_3 (gsl_matrix *G, gsl_matrix *G_inv);
void block_inv_4 (gsl_matrix *G, gsl_matrix *G_inv);

void time_elapsed (time_t &a,time_t &b);

inline double _max (double a, double b) { return (a > b ? a : b); };
inline double _min (double a, double b) { return (a > b ? a : b); };
inline int _abs (int a) { return (a > 0 ? a : -a); }; 

void usage() {
  printf("\t\t\t\t\t\t\t" SIGNATURE "\n");
  printf("  BE2 (Boundary Element method for 2 surfaces)\n");
  printf("\n");
  printf("  calculate translational/rotational diffusion tensors\n");
  printf("\n");
  printf("  Usage: be2 [-msms|-gts|-off] file_prefix [options]\n");
  printf("\n");
  printf("  options:\n");
  printf("    [-rR gts_prefix]  read reference surface  (GTS format)\n");
  printf("    [-t #]            temperature (default: 293.15 K)\n");
  printf("    [-v #]            viscosity (default: 1.002 cP)\n");
  printf("    [-g #]            generous boundary (default: 6 A)\n");
  printf("    [-e #]            epsilon (default: 22 A)\n");
  printf("    [-gmap # # #]     map of generous boundary (min max step)\n");
  printf("    [-emap # # #]     map of epsilon (min max step)\n");
  printf("    [-coff]           OFF colormap output of velocity factor\n");
  printf("    [-a 1|2|3]        algorithm for velocity factor calc. (default: 1)\n");
  printf("    [-rP pdb_prefix]  (algorithm 2|3) read sample PDB coordinates\n");
  printf("    [-r # # #]        (algorithm 2|3) rigid segment chainId, N-, C-ter.\n");
  printf("    [-verb]           long output\n");
  printf("    [-h]              display this message\n");
  printf("\n");
  exit (-1);
  }

/*
  N		  : number of surface elements
  dSVert	  : vertice coord. of surface elements (N x 9)
  dSSide	  : side vectors of surface elements (N x 9)
  dSLength	  : side vector lengths of surface elements (N x 3)
  dS		  : area of surface elements (N x 1)
  Ci		  : incenter (N x 3)
  C		  : centroid (N x 3)
  G		  : Geometry matrix (3N x 3N)
  G_inv		  : inverse matrix of G (3N x 3N)
  OR		  : center of diffusion (3 x 1)
  vf		  : velocity factor (N x 1)
  NR		  : number of reference surface elements
  RC		  : centroid of reference surface elements (NR x 3)
  dRVert	  : vertices coord. of reference surface elements (NR x 9)
  Ktt,Ktr,Krt,Krr : friction tensors (3 x 3)
  Dtt,Dtr,Drt,Drr : diffusion tensors (3 x 3)
*/

int main(int argc, char *argv[]) 
{
  /* checking elapsed time */
  time_t start,stage,finish;
  time(&start);
  
  extern double temperature, viscosity, epsilon, gboundary;
  extern bool verb;

  size_t N,NR,M,i,j,k,u,v,delta; 
  int signum,status;
  double A[3],l[3],h[3],SI[3],CI[3];
  double a,b,c,r0,sc,cc,dp,cosA,ln;
  double d1,d2,d3;
  double abserr,T0,S0,S1,S2,S12,S11,S22;

  char  dS_file_prefix[MAXLENGTH];
  char  dR_file_prefix[MAXLENGTH];
  char pdb_file_prefix[MAXLENGTH];

  gsl_matrix *dSVert = NULL;
  gsl_matrix *dRVert = NULL;
  gsl_matrix *RC = NULL;  

  Attrib dSa;
  vector < Atom > coor;

  N = 0;
  size_t dS_format = _GTS_FORMAT_;
  bool dR_flag = false;
  bool coor_flag = false;

  if (argc < 2) usage();

  i = 1;
  while (i < (size_t)argc) {
    if (strcmp(argv[i],"-help") == 0 || strcmp(argv[i],"-h") == 0) 
      usage();
    if (strcmp(argv[i],"-gts" ) == 0) {
      dS_format = _GTS_FORMAT_;
      strcpy(dS_file_prefix,argv[++i]);
      }
    if (strcmp(argv[i],"-msms") == 0) {
      dS_format = _MSMS_FORMAT_;
      strcpy(dS_file_prefix,argv[++i]);
      }
    if (strcmp(argv[i],"-off") == 0) {
      dS_format = _OFF_FORMAT_;
      strcpy(dS_file_prefix,argv[++i]);
      }
    if (strcmp(argv[i],"-rR") == 0) {
      dR_flag = true;
      strcpy(dR_file_prefix,argv[++i]);
      }
    if (strcmp(argv[i],"-rP") == 0) {
      coor_flag = true;
      strcpy(pdb_file_prefix,argv[++i]);
      } // -rP

    if (strcmp(argv[i],"-r") == 0) {
      sscanf(argv[++i],"%c" ,&Rchain_);
      sscanf(argv[++i],"%d" ,&RNter_);
      sscanf(argv[++i],"%d" ,&RCter_);
      Rchain.push_back (Rchain_);
      RNter.push_back (RNter_);
      RCter.push_back (RCter_);
      } // -r

    if (strcmp(argv[i],"-a") == 0) sscanf(argv[++i],"%d" ,&algorithm);
    if (strcmp(argv[i],"-t") == 0) sscanf(argv[++i],"%lf",&temperature);
    if (strcmp(argv[i],"-v") == 0) sscanf(argv[++i],"%lf",&viscosity);
    if (strcmp(argv[i],"-e") == 0) sscanf(argv[++i],"%lf",&epsilon);
    if (strcmp(argv[i],"-g") == 0) sscanf(argv[++i],"%lf",&gboundary);
    if (strcmp(argv[i],"-gmap") == 0) {
      gmap = true;
      sscanf(argv[++i],"%lf",&map_g_min);
      sscanf(argv[++i],"%lf",&map_g_max);
      sscanf(argv[++i],"%lf",&map_g_step);
      }
    if (strcmp(argv[i],"-emap") == 0) {
      emap = true;
      sscanf(argv[++i],"%lf",&map_e_min);
      sscanf(argv[++i],"%lf",&map_e_max);
      sscanf(argv[++i],"%lf",&map_e_step);
      }
    if (strcmp(argv[i],"-coff") == 0) coff = true;
    if (strcmp(argv[i],"-verb") == 0) verb = true;
    i++;
    }

  time(&stage);

  /* constants */
  printf("Temperature  : %g K\n",temperature);
  printf("Viscosity    : %g cP\n",viscosity);
  if (dR_flag)
  printf("Algorithm    : %d\n",algorithm);

  /* read sample surface elements */
  if (dS_format == _GTS_FORMAT_ ) {
    N = read_gts_surface (dS_file_prefix,NULL);
    dSVert = gsl_matrix_alloc (N, 9);
    read_gts_surface (dS_file_prefix,dSVert);
    calc_dimension (dSVert);
    }
  if (dS_format == _MSMS_FORMAT_ ) {
    N = read_msms_surface (dS_file_prefix,NULL);
    dSVert = gsl_matrix_alloc (N, 9);
    read_msms_surface (dS_file_prefix,dSVert);
    calc_dimension (dSVert);
    }
  if (dS_format == _OFF_FORMAT_ ) {
    N = read_geom_surface(dS_file_prefix,NULL);
    dSVert = gsl_matrix_alloc (N, 9);
    read_geom_surface (dS_file_prefix,dSVert);
    calc_dimension (dSVert);
    }

  if (N == 0) usage();

  /* read reference surface elements */
  if (dR_flag) {
    NR = read_gts_surface (dR_file_prefix,NULL);
    dRVert = gsl_matrix_alloc (NR, 9);
    RC = gsl_matrix_alloc (NR,3);
    read_gts_surface (dR_file_prefix,dRVert);
    calc_dimension (dRVert);
    calc_surf_geometry (dRVert,RC,NULL,NULL,NULL);
    }

  /* read pdb coordiantes */
  if (coor_flag) {
    M = read_pdb_coor (pdb_file_prefix, coor);
    printf ("PDB coord.   : %s.pdb (Atoms: %d)\n",pdb_file_prefix,M);
    if (0) {
      for (j=0;j<M;j++) {
	printf ("%5d %3s %4d %4s %8.3f %8.3f %8.3f\n",
	  (j+1), coor[j].resn.c_str(),
	  coor[j].r, coor[j].name.c_str(),
	  coor[j].x, coor[j].y, coor[j].z);
	} // coordinate
      } // debug-only
    }

  scalar_product D;
  gsl_function F;
  gsl_matrix_view Gkj; 
  gsl_vector_view v1,v2,v3,va,vb,vc,ck,vj,sj,gv;
  gsl_permutation *p  = gsl_permutation_alloc (3);
  gsl_matrix *G       = gsl_matrix_alloc (3*N,3*N);
  gsl_vector *dS      = gsl_vector_alloc (N);
  gsl_matrix *dSSide  = gsl_matrix_alloc (N,9);
  gsl_matrix *dSLength= gsl_matrix_alloc (N,3);
  gsl_matrix *C       = gsl_matrix_alloc (N,3);
  gsl_matrix *G_lab   = gsl_matrix_alloc (3,3);
  gsl_matrix *G_loc   = gsl_matrix_alloc (3,3);
  gsl_matrix *R       = gsl_matrix_alloc (3,3);
  gsl_matrix *R_LU    = gsl_matrix_alloc (3,3);
  gsl_matrix *R_inv   = gsl_matrix_alloc (3,3);
  gsl_matrix *temp    = gsl_matrix_alloc (3,3);
  gsl_vector *Ci      = gsl_vector_alloc (3);
  gsl_vector *n       = gsl_vector_alloc (3);
  gsl_vector *n_      = gsl_vector_alloc (3);
  gsl_vector *va_     = gsl_vector_alloc (3);
  gsl_vector *vb_     = gsl_vector_alloc (3);
  gsl_vector *vc_     = gsl_vector_alloc (3);
  gsl_vector *y_      = gsl_vector_alloc (3);
  gsl_vector *vt      = gsl_vector_alloc (3);
  gsl_vector *vd1     = gsl_vector_alloc (3);
  gsl_vector *vd2     = gsl_vector_alloc (3);
  gsl_vector *vd3     = gsl_vector_alloc (3);
  gsl_vector *vf      = gsl_vector_alloc (N);


/* Calculate surface geometries: centroid, side vectors, length, area */

  calc_surf_geometry (dSVert,C,dSSide,dSLength,dS);

/* Calculate surface element attributes */
  if (RC != NULL)
    calc_attrib (C, RC, coor, dSa);

/*
  Calculate 3N x 3N G matrix
  G matrix contains geometric information relevant to the
  hydrodynamic interaction calculation
  Place the coordinate xj at the center of the small patch Dj and 
  take the surface stress force f(x) to be a constant over the 
  entire patch area. 
  Gkj = integration T(x,yk)dSx
  Hydrodynamic interaction felt at k surface by j surface.
*/

  gsl_integration_workspace *wqag = 
	gsl_integration_workspace_alloc (QAG_PREC_INTVL);

  for (k=0; k<N; k++) {
    for (j=0; j<N; j++) {

      /* common for diagonal and off-diagonal calculations */
      /* vertices, j */
      vj = gsl_matrix_row (dSVert,j);
      v1 = gsl_vector_subvector (&vj.vector,0,3);
      v2 = gsl_vector_subvector (&vj.vector,3,3);
      v3 = gsl_vector_subvector (&vj.vector,6,3);
      /* sides, j */
      sj = gsl_matrix_row (dSSide,j);
      va = gsl_vector_subvector (&sj.vector,0,3);
      vb = gsl_vector_subvector (&sj.vector,3,3);
      vc = gsl_vector_subvector (&sj.vector,6,3);
      /* lengths, j */
      a = gsl_matrix_get (dSLength,j,0);
      b = gsl_matrix_get (dSLength,j,1);
      c = gsl_matrix_get (dSLength,j,2);

      if (j==k) { /* DIAGONAL ELEMENTS */
	/* unit vectors */
	gsl_vector_set_zero(va_);
	gsl_vector_set_zero(vb_);
	gsl_vector_set_zero(vc_);
	gsl_vector_set_zero(n_);
	gsl_blas_daxpy (1/a,&va.vector,va_);
	gsl_blas_daxpy (1/b,&vb.vector,vb_);
	gsl_blas_daxpy (1/c,&vc.vector,vc_);
	gsl_blas_ddot(vb_,vc_,&cosA); 
	A[0] = acos( cosA); /* angle A */
	gsl_blas_ddot(va_,vc_,&cosA);
	A[1] = acos(-cosA); /* angle B */
	gsl_blas_ddot(va_,vb_,&cosA); 
	A[2] = acos( cosA); /* angle C */
	/* n = vb x va */
	gsl_vector_set(n,0,
	  gsl_vector_get(&vb.vector,1)*gsl_vector_get(&va.vector,2)-
	  gsl_vector_get(&vb.vector,2)*gsl_vector_get(&va.vector,1));
	gsl_vector_set(n,1,
	  gsl_vector_get(&vb.vector,2)*gsl_vector_get(&va.vector,0)-
	  gsl_vector_get(&vb.vector,0)*gsl_vector_get(&va.vector,2));
	gsl_vector_set(n,2,
	  gsl_vector_get(&vb.vector,0)*gsl_vector_get(&va.vector,1)-
	  gsl_vector_get(&vb.vector,1)*gsl_vector_get(&va.vector,0));
	/* ln = 2*Area */
	ln = gsl_blas_dnrm2(n); 
	gsl_blas_daxpy (1/ln,n,n_);
	/* radius of incircle */
	r0 = ln/(a+b+c);
	/* incenter, Ci */
	gsl_vector_set_zero(Ci);
	gsl_blas_daxpy (a/(a+b+c),&v3.vector,Ci);
	gsl_blas_daxpy (b/(a+b+c),&v2.vector,Ci);
	gsl_blas_daxpy (c/(a+b+c),&v1.vector,Ci);
	/* distance from incenter to vertices */
	gsl_vector_memcpy(vt,&v3.vector);
	gsl_vector_sub(vt,Ci);	
	l[0] = gsl_blas_dnrm2(vt); /* la */
	gsl_vector_memcpy(vt,&v2.vector);
	gsl_vector_sub(vt,Ci);	
	l[1] = gsl_blas_dnrm2(vt); /* lb */
	gsl_vector_memcpy(vt,&v1.vector);
	gsl_vector_sub(vt,Ci);	
	l[2] = gsl_blas_dnrm2(vt); /* lc */
	/* y-axis, y_ */
	gsl_vector_set(y_,0,
	  gsl_vector_get(n_,1)*gsl_vector_get(va_,2)-
	  gsl_vector_get(n_,2)*gsl_vector_get(va_,1));
	gsl_vector_set(y_,1,
	  gsl_vector_get(n_,2)*gsl_vector_get(va_,0)-
	  gsl_vector_get(n_,0)*gsl_vector_get(va_,2));
	gsl_vector_set(y_,2,
	  gsl_vector_get(n_,0)*gsl_vector_get(va_,1)-
	  gsl_vector_get(n_,1)*gsl_vector_get(va_,0));
	/* rotation matrix */
	gv = gsl_matrix_row (R,0);
	gsl_vector_memcpy(&gv.vector,va_);
	gv = gsl_matrix_row (R,1);
	gsl_vector_memcpy(&gv.vector,y_);
	gv = gsl_matrix_row (R,2);
	gsl_vector_memcpy(&gv.vector,n_);
	/* inverse matrix of R */
	gsl_permutation_init (p);
	gsl_matrix_memcpy(R_LU,R);
	gsl_linalg_LU_decomp (R_LU,p,&signum);
	gsl_linalg_LU_invert (R_LU,p,R_inv);
	/* for integration of outside of incircle */
	dp = 0.0;
	for (i=0;i<3;i++) {
	  h[i]  = l[i]/r0;
	  CI[i] = 3.0 -h[i] -2.0/h[i];
	  SI[i] = 2.0*(log(h[i]+sqrt(h[i]*h[i]-1))-sqrt(1-1/(h[i]*h[i])));
	  dp += (M_PI-A[i])*(h[i]-1)-2.0*(h[i]*acos(1/h[i])
	    -log(h[i]+sqrt(h[i]*h[i]-1.0)));
	  }
	sc =-0.5*cos(A[1]-A[2])*(SI[0]*cos(A[0])+CI[0]*sin(A[0]));
	sc+= 0.25*(SI[1]*(1+cos(2*A[1]))+CI[1]*sin(2*A[1]));
	sc+= 0.25*(SI[2]*(1+cos(2*A[2]))+CI[2]*sin(2*A[2]));
	cc =-(sin(A[1]-A[2])*cos(A[0])*CI[0]+sin(A[1]-A[2])*sin(A[0])*SI[0]);
	cc+= (sin(A[1])*sin(A[1])*CI[1]+sin(A[1])*cos(A[1])*SI[1]);
	cc+=-(sin(A[2])*sin(A[2])*CI[2]+sin(A[2])*cos(A[2])*SI[2]);

	/* integration in local coordinate */
	gsl_matrix_set(G_loc,0,0,3.0*M_PI+1.5*dp+0.25*sc);
	gsl_matrix_set(G_loc,0,1,0.5*cc);
	gsl_matrix_set(G_loc,0,2,0);
	gsl_matrix_set(G_loc,1,0,0.5*cc);
	gsl_matrix_set(G_loc,1,1,3.0*M_PI+1.5*dp-0.25*sc);
	gsl_matrix_set(G_loc,1,2,0);
	gsl_matrix_set(G_loc,2,0,0);
	gsl_matrix_set(G_loc,2,1,0);
	gsl_matrix_set(G_loc,2,2,2.0*M_PI+dp);
	gsl_matrix_scale(G_loc,r0);
	/* transform to laboratory coordinate (R_inv G_loc) R ,(AB)C = A(BC) */
	gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,R_inv,G_loc,0,temp);
	gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,temp,R,0,G_lab);
	} // if j==k

  
      /* OFF-DIAGONAL ELEMENTS */
      if (j!=k) { 
	/* 
	  Cc: centroid of k surface
	  v1,v2,v3: vertices of j surface
	  va  = v2 - v1
	  vb  = v3 - v1
	  vc  = v3 - v2
	  vd1 = v1 - Cc
	  vd2 = v2 - Cc
	  vd3 = v3 - Cc
	*/
	ck = gsl_matrix_row (C,k); /* k, centroid */
	gsl_vector_memcpy(vd1,&v1.vector);
	gsl_vector_sub(vd1,&ck.vector); 
	d1 = gsl_blas_dnrm2(vd1);
	gsl_vector_memcpy(vd2,&v2.vector);
	gsl_vector_sub(vd2,&ck.vector);
	d2 = gsl_blas_dnrm2(vd2);
	gsl_vector_memcpy(vd3,&v3.vector);
	gsl_vector_sub(vd3,&ck.vector);
	d3 = gsl_blas_dnrm2(vd3);

	/* scalar product for integration */
	gsl_blas_ddot (&va.vector,vd1,&D.ad1);
	gsl_blas_ddot (&va.vector,vd2,&D.ad2);
	gsl_blas_ddot (&va.vector,&vc.vector,&D.ac);
	gsl_blas_ddot (&va.vector,&vb.vector,&D.ab);
	gsl_blas_ddot (&vb.vector,&vc.vector,&D.bc);
	gsl_blas_ddot (&vb.vector,vd1,&D.bd1);
	gsl_blas_ddot (&vb.vector,vd3,&D.bd3);
	gsl_blas_ddot (&vc.vector,vd2,&D.cd2);
	gsl_blas_ddot (&vc.vector,vd3,&D.cd3);
	gsl_blas_ddot (vd1,vd2,&D.d1d2);
	D.a = a;
	D.sqa = a*a;
	D.sqb = b*b;
	D.sqc = c*c;
	D.sqd1 = d1*d1;
	D.sqd2 = d2*d2;
	D.sqd3 = d3*d3;

	/* numerical integration by QAG */
	gsl_set_error_handler_off();
	F.params = &D;
	F.function = &dyt0;
	status = gsl_integration_qag (&F,0,1,EPSABS,EPSREL,
		  QAG_PREC_INTVL,QAG_RULE,wqag,&T0,&abserr);
	if (status) integ_err(k,j,D);
	F.function = &dyg0;
	status = gsl_integration_qag (&F,0,1,EPSABS,EPSREL,
		QAG_PREC_INTVL,QAG_RULE,wqag,&S0,&abserr);
	if (status) integ_err(k,j,D);
	F.function = &dyg1;
	status = gsl_integration_qag (&F,0,1,EPSABS,EPSREL,
		QAG_PREC_INTVL,QAG_RULE,wqag,&S1,&abserr);
	if (status) integ_err(k,j,D);
	F.function = &dyyg0;
	status = gsl_integration_qag (&F,0,1,EPSABS,EPSREL,
		QAG_PREC_INTVL,QAG_RULE,wqag,&S2,&abserr);
	if (status) integ_err(k,j,D);
	F.function = &dyyg1;
	status = gsl_integration_qag (&F,0,1,EPSABS,EPSREL,
		QAG_PREC_INTVL,QAG_RULE,wqag,&S12,&abserr);
	if (status) integ_err(k,j,D);
	F.function = &dxxxg2;
	status = gsl_integration_qag (&F,0,1,EPSABS,EPSREL,
		QAG_PREC_INTVL,QAG_RULE,wqag,&S11,&abserr);
	if (status) integ_err(k,j,D);
	F.function = &dyyyg0;
	status = gsl_integration_qag (&F,0,1,EPSABS,EPSREL,
		QAG_PREC_INTVL,QAG_RULE,wqag,&S22,&abserr);
	if (status) integ_err(k,j,D);
	gsl_set_error_handler(NULL);

	/* G_lab */
	for(u=0; u<3; u++) {
	  for (v=0; v<3; v++) {
	    if (u==v) delta = 1; else delta = 0; 
	    gsl_matrix_set(G_lab,u,v,
	      2*gsl_vector_get(dS,j)*(
	      T0*delta + S0*gsl_vector_get(vd1,u)*gsl_vector_get(vd1,v) +
	      S1*(gsl_vector_get(vd1,u)*gsl_vector_get(&va.vector,v)+
		  gsl_vector_get(&va.vector,u)*gsl_vector_get(vd1,v)) +
	      S2*(gsl_vector_get(vd1,u)*gsl_vector_get(&vb.vector,v)+
		  gsl_vector_get(&vb.vector,u)*gsl_vector_get(vd1,v)) +
	      S12*(gsl_vector_get(&va.vector,u)*gsl_vector_get(&vb.vector,v)+
		  gsl_vector_get(&vb.vector,u)*gsl_vector_get(&va.vector,v))+
	      S11*gsl_vector_get(&va.vector,u)*gsl_vector_get(&va.vector,v) +
	      S22*gsl_vector_get(&vb.vector,u)*gsl_vector_get(&vb.vector,v)));
	    } // v
	  } // u

	} // if j!=k
      /* fill-in 3N x 3N supermatrix G */
      Gkj = gsl_matrix_submatrix (G, 3*k, 3*j, 3, 3);
      gsl_matrix_memcpy(&Gkj.matrix, G_lab);
      } // j loop
    } // k loop
  gsl_matrix_free (G_lab);
  gsl_matrix_free (G_loc);
  gsl_matrix_free (R);
  gsl_matrix_free (R_LU);
  gsl_matrix_free (R_inv);
  gsl_vector_free (Ci);
  gsl_vector_free (va_);
  gsl_vector_free (vb_);
  gsl_vector_free (vc_);
  gsl_vector_free (n);
  gsl_vector_free (n_);
  gsl_vector_free (y_);
  gsl_vector_free (vt);
  gsl_vector_free (vd1);
  gsl_vector_free (vd2);
  gsl_vector_free (vd3);
  time(&finish);
  printf("Calculating G matrix (%4d x %4d) ...... ",3*N,3*N);
  time_elapsed (stage, finish);

/* calculate inverse matrix of G using LU decomposition */

  time (&stage);

  gsl_matrix *G_inv   = gsl_matrix_alloc (G->size1,G->size2);

/* quadruple block inversion */

  block_inv (G,G_inv); 

/*
  gsl_permutation *q  = gsl_permutation_alloc (G->size1);
  gsl_permutation_init (q);
  gsl_linalg_LU_decomp (G,q,&signum);
  gsl_linalg_LU_invert (G,q,G_inv);
  gsl_permutation_free (q);
*/

  gsl_matrix_free (G);

  time (&finish);
  printf("Inverting   G matrix (%4d x %4d) ...... ",3*N,3*N);
  time_elapsed (stage, finish);

/*
  surface area statistics
*/
  double dSmin,dSmax,Area;
  for (Area=0,k=0;k<N;k++)
    Area += gsl_vector_get(dS,k);
  gsl_vector_minmax(dS,&dSmin,&dSmax);
  printf("Surface_Area (%10.6f ... %10.6f) Sum= %g A^2\n",dSmin,dSmax,Area);

/*
  calculate friction tensor & diffusion Tensor
  update center of diffusion & recalculate 
  note: Dtt and Dtr are affected by the center of diffusion
*/

  time (&stage);

  gsl_matrix *Ktt     = gsl_matrix_alloc (3,3);
  gsl_matrix *Ktr     = gsl_matrix_alloc (3,3);
  gsl_matrix *Krt     = gsl_matrix_alloc (3,3);
  gsl_matrix *Krr     = gsl_matrix_alloc (3,3);
  gsl_matrix *Drr     = gsl_matrix_alloc (3,3);
  gsl_matrix *Dtr     = gsl_matrix_alloc (3,3);
  gsl_matrix *Drt     = gsl_matrix_alloc (3,3);
  gsl_matrix *Dtt     = gsl_matrix_alloc (3,3);
  gsl_vector *OR      = gsl_vector_alloc (3);
  gsl_matrix *evec    = gsl_matrix_alloc (3,3);
  gsl_vector *eval    = gsl_vector_alloc (3);
  gsl_eigen_symmv_workspace * w 
		      = gsl_eigen_symmv_alloc (3);

  if (!gmap) map_g_min = map_g_max = gboundary;
  if (!emap) map_e_min = map_e_max = epsilon;
  
  gboundary = map_g_min;
  while (gboundary <= map_g_max) {
  epsilon = map_e_min;
  while (epsilon <= map_e_max) {

  /* set initial velocity factor to 1 */
  gsl_vector_set_all (vf, 1.0);

  /* velocity field assumption for soft surface */
  /* calculate velocity factor if reference surface is given */
  if (RC != NULL) {
    calc_velocity_factor (C,dSa,vf);
    if (coff) write_geom_surface (dS_file_prefix,dSVert,vf,0); /* colormap */
    }
    
  /* set initial point to (0,0,0) */
  gsl_vector_set_zero (OR);
  calc_friction_tensor (G_inv,C,OR,dS,vf,Ktt,Ktr,Krt,Krr);
  calc_diffusion_tensor(Ktt,Ktr,Krt,Krr,Dtt,Dtr,Drt,Drr);
  calc_center_of_diffusion (Drr,Dtr,OR);

  printf("\nCenter_of_Diffusion %.4f %.4f %.4f\n",
    -gsl_vector_get(OR,0), 
    -gsl_vector_get(OR,1),
    -gsl_vector_get(OR,2));

  /* recalculate friction and diffusion tensor with adjusted OR */
  calc_friction_tensor (G_inv,C,OR,dS,vf,Ktt,Ktr,Krt,Krr);
  calc_diffusion_tensor(Ktt,Ktr,Krt,Krr,Dtt,Dtr,Drt,Drr);

  if (verb)
    display_matrix (Ktt,Ktr,Krt,Krr,Dtt,Dtr,Drt,Drr);

/*
  Eigenvalues and Eigenvectors
*/
  gsl_eigen_symmv (Dtt,eval,evec,w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  printf("Dtt Eigenvalues (10^-7 cm^2 s^-1) | Eigenvectors\n");
  for (u=0;u<3;u++) {
    printf("Dtt %1d %10.3f | ",u+1,
      C_TRANSLATION*temperature/viscosity*gsl_vector_get(eval,u));
    gsl_vector_view evec_u = gsl_matrix_column (evec,u);
    printf("%10.3e %10.3e %10.3e\n",
      gsl_matrix_get(evec,0,u),
      gsl_matrix_get(evec,1,u),
      gsl_matrix_get(evec,2,u));
    }
  gsl_eigen_symmv (Drr,eval,evec,w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  printf("Drr Eigenvalues (10^7 s^-1) | Eigenvectors\n");
  double Diso = 0.0;
  for (u=0;u<3;u++) {
    printf("Drr %1d %10.3f | ",u+1,
      C_ROTATION*temperature/viscosity*gsl_vector_get(eval,u));
    printf("%10.3e %10.3e %10.3e\n",
      gsl_matrix_get(evec,0,u),
      gsl_matrix_get(evec,1,u),
      gsl_matrix_get(evec,2,u));
    Diso += C_ROTATION*temperature/viscosity*gsl_vector_get(eval,u);
    }
  if (RC!=NULL) 
    printf("# gamma %g eps %g 1/(6Diso) %7.3f ns\n\n",
	gboundary,epsilon,50.0/Diso);
  else
    printf("# rigid 1/(6Diso) %7.3f ns\n\n",50.0/Diso);

  epsilon += map_e_step;
  } /* loop eps */
  gboundary += map_g_step;
  } /* loop generous boundary */

  /* 
  supplementary results based on rigid body
  */

  if (RC != NULL && 0) {

    /* reset velocity factor */
    gsl_vector_set_all (vf, 1.0);

    /* set initial point to (0,0,0) */
    gsl_vector_set_zero (OR);
    calc_friction_tensor (G_inv,C,OR,dS,vf,Ktt,Ktr,Krt,Krr);
    calc_diffusion_tensor(Ktt,Ktr,Krt,Krr,Dtt,Dtr,Drt,Drr);
    calc_center_of_diffusion (Drr,Dtr,OR);
    printf("\nCenter_of_Diffusion %.4f %.4f %.4f\n",
      -gsl_vector_get(OR,0), 
      -gsl_vector_get(OR,1), 
      -gsl_vector_get(OR,2));

    /* recalculate friction and diffusion tensor with updated OR */
    calc_friction_tensor (G_inv,C,OR,dS,vf,Ktt,Ktr,Krt,Krr);
    calc_diffusion_tensor(Ktt,Ktr,Krt,Krr,Dtt,Dtr,Drt,Drr);
    if (verb) 
	display_matrix (Ktt,Ktr,Krt,Krr,Dtt,Dtr,Drt,Drr);

    /* Eigenvalues and Eigenvectors */
    gsl_eigen_symmv (Dtt,eval,evec,w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    printf("Dtt Eigenvalues (10^-7 cm^2 s^-1) | Eigenvectors\n");
    for (u=0;u<3;u++) {
      printf("Dtt %1d %10.3f | ",u+1,
	C_TRANSLATION*temperature/viscosity*gsl_vector_get(eval,u));
      gsl_vector_view evec_u = gsl_matrix_column (evec,u);
      printf("%10.3e %10.3e %10.3e\n",
	gsl_matrix_get(evec,0,u),
	gsl_matrix_get(evec,1,u),
	gsl_matrix_get(evec,2,u));
      }
    gsl_eigen_symmv (Drr,eval,evec,w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    printf("Drr Eigenvalues (10^7 s^-1) | Eigenvectors\n");
    double Diso = 0.0;
    for (u=0;u<3;u++) {
      printf("Drr %1d %10.3f | ",u+1,
	C_ROTATION*temperature/viscosity*gsl_vector_get(eval,u));
      printf("%10.3e %10.3e %10.3e\n",
	gsl_matrix_get(evec,0,u),
	gsl_matrix_get(evec,1,u),
	gsl_matrix_get(evec,2,u));
      Diso += C_ROTATION*temperature/viscosity*gsl_vector_get(eval,u);
      }
    printf("# rigid 1/(6Diso) %7.3f ns\n\n",50.0/Diso);
    } /* if RC is defined */

  gsl_eigen_symmv_free (w);
  gsl_vector_free (eval);
  gsl_matrix_free (evec);

  /* freeing memory */
  gsl_matrix_free (G_inv);
  gsl_vector_free (dS);
  gsl_matrix_free (dSVert);
  gsl_matrix_free (C);
  gsl_matrix_free (temp);
  gsl_matrix_free (Ktt);
  gsl_matrix_free (Ktr);
  gsl_matrix_free (Krt);
  gsl_matrix_free (Krr);
  gsl_matrix_free (Dtt);
  gsl_matrix_free (Dtr);
  gsl_matrix_free (Drt);
  gsl_matrix_free (Drr);
  gsl_vector_free (OR);
  gsl_vector_free (vf);
  gsl_permutation_free (p);

/*
  time(&finish);
  printf("Calculating Tensors ...... ");
  time_elapsed (stage, finish);
*/

  /* Total elapsed time */
  if (verb) {
    time(&finish);
    printf("\n");
    printf("started  %s",ctime(&start));
    printf("finished %s",ctime(&finish));
    printf("Total ");
    time_elapsed (start, finish);
    }
}// main()


void display_matrix (
    gsl_matrix *Ktt, gsl_matrix *Ktr,
    gsl_matrix *Krt, gsl_matrix *Krr,
    gsl_matrix *Dtt, gsl_matrix *Dtr,
    gsl_matrix *Drt, gsl_matrix *Drr)
{
    int u;
    printf("\nKtt Ktr :\n");
    for (u=0;u<3;u++) {
	printf("| %10.3e %10.3e %10.3e |   | %10.3e %10.3e %10.3e |\n",
	    gsl_matrix_get(Ktt,u,0), 
	    gsl_matrix_get(Ktt,u,1), 
	    gsl_matrix_get(Ktt,u,2), 
	  gsl_matrix_get(Ktr,u,0),
	  gsl_matrix_get(Ktr,u,1),
	  gsl_matrix_get(Ktr,u,2));
	    }
	printf("Krt Krr :\n");
      for (u=0;u<3;u++) {
	printf("| %10.3e %10.3e %10.3e |   | %10.3e %10.3e %10.3e |\n",
	  gsl_matrix_get(Krt,u,0), 
	  gsl_matrix_get(Krt,u,1), 
	  gsl_matrix_get(Krt,u,2), 
	  gsl_matrix_get(Krr,u,0),
	  gsl_matrix_get(Krr,u,1),
	  gsl_matrix_get(Krr,u,2));
	}
      printf("\nDtt Dtr :\n");
      for (u=0;u<3;u++) {
	printf("| %10.3e %10.3e %10.3e |   | %10.3e %10.3e %10.3e |\n",
	  gsl_matrix_get(Dtt,u,0),
	  gsl_matrix_get(Dtt,u,1),
	  gsl_matrix_get(Dtt,u,2),
	  gsl_matrix_get(Dtr,u,0),
	  gsl_matrix_get(Dtr,u,1),
	  gsl_matrix_get(Dtr,u,2));
	}
      printf("Drt Drr :\n");
      for (u=0;u<3;u++) {
	printf("| %10.3e %10.3e %10.3e |   | %10.3e %10.3e %10.3e |\n",
	  gsl_matrix_get(Drt,u,0),
	  gsl_matrix_get(Drt,u,1),
	  gsl_matrix_get(Drt,u,2),
	  gsl_matrix_get(Drr,u,0),
	  gsl_matrix_get(Drr,u,1),
	  gsl_matrix_get(Drr,u,2));
	}
}


/*
  calculate friction tensor (K) using force and torque relations

  (1) for flexible surface patches, 
      velocity factor = d*exp(-d/EPS)
      in which d is the distance between sample surface patch and 
      reference surface patch to account for the velocity and
      angular velocity decay from the reference surface.

  (2) cross product as matrix multiplication 
      general property of cross product: a x b = - b x a
      a (x,y,z) x  (...) =  A (...)
      (...) x a (x,y,z)  = -a (x,y,z) x (...) = -A (...)
      A = | 0 -z  y |
	  | z  0 -x |
	  |-y  x  0 |
	
  (3) gsl_blas_dgemm: (opA,opB,alpha,A,B,beta,C)
      sum C = \alpha op(A) op(B) + \beta C
*/

void calc_friction_tensor (
  gsl_matrix *G_inv,
  gsl_matrix *C,
  gsl_vector *OR,
  gsl_vector *dS,
  gsl_vector *vf,
  gsl_matrix *Ktt,
  gsl_matrix *Ktr,
  gsl_matrix *Krt,
  gsl_matrix *Krr) 
{
  size_t j,k,N = dS->size;
  double x,y,z,vfk;
  gsl_matrix *rjX  = gsl_matrix_alloc (3,3);
  gsl_matrix *rkX  = gsl_matrix_alloc (3,3);
  gsl_matrix *E    = gsl_matrix_alloc (3,3);
  gsl_matrix *temp = gsl_matrix_alloc (3,3);
  gsl_matrix_view iGjk;

  gsl_matrix_set_identity (E);
  gsl_matrix_set_zero (Ktt);
  gsl_matrix_set_zero (Ktr);
  gsl_matrix_set_zero (Krt);
  gsl_matrix_set_zero (Krr);
  
  for (j=0;j<N;j++) {
    x  = gsl_matrix_get(C,j,0) + gsl_vector_get(OR,0);
    y  = gsl_matrix_get(C,j,1) + gsl_vector_get(OR,1);
    z  = gsl_matrix_get(C,j,2) + gsl_vector_get(OR,2);
    /* rj x , cross product as matrix multiplication */
    gsl_matrix_set (rjX,0,0, 0);
    gsl_matrix_set (rjX,0,1,-z);
    gsl_matrix_set (rjX,0,2, y);
    gsl_matrix_set (rjX,1,0, z);
    gsl_matrix_set (rjX,1,1, 0);
    gsl_matrix_set (rjX,1,2,-x);
    gsl_matrix_set (rjX,2,0,-y);
    gsl_matrix_set (rjX,2,1, x);
    gsl_matrix_set (rjX,2,2, 0);
    for (k=0;k<N;k++) {
      x  = gsl_matrix_get(C,k,0) + gsl_vector_get(OR,0);
      y  = gsl_matrix_get(C,k,1) + gsl_vector_get(OR,1);
      z  = gsl_matrix_get(C,k,2) + gsl_vector_get(OR,2);
      /* rk x , cross product as matrix multiplication */
      gsl_matrix_set (rkX,0,0, 0);
      gsl_matrix_set (rkX,0,1,-z);
      gsl_matrix_set (rkX,0,2, y);
      gsl_matrix_set (rkX,1,0, z);
      gsl_matrix_set (rkX,1,1, 0);
      gsl_matrix_set (rkX,1,2,-x);
      gsl_matrix_set (rkX,2,0,-y);
      gsl_matrix_set (rkX,2,1, x);
      gsl_matrix_set (rkX,2,2, 0);

      iGjk = gsl_matrix_submatrix(G_inv, 3*j, 3*k, 3, 3);
      vfk = gsl_vector_get(vf,k); /* velocity factor */

      /* Ktt += dSj*iGjk */
      gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,
	    vfk*gsl_vector_get(dS,j),&iGjk.matrix,E,1,Ktt);

      /* Ktr += -dSj*iGjk rkX */
      gsl_blas_dgemm (CblasNoTrans,CblasNoTrans, 
	    -vfk*gsl_vector_get(dS,j),&iGjk.matrix,rkX,1,Ktr);

      /* Krt += dSj*(rjX iGjk) */
      gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,
	    vfk*gsl_vector_get(dS,j),rjX,&iGjk.matrix,1,Krt);

      /* Krr += -dSj*(rjX iGjk rkX) */
      gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,
	    -vfk*gsl_vector_get(dS,j),rjX,&iGjk.matrix,0,temp);
      gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,temp,rkX,1,Krr);
      } // k
    } // j
  gsl_matrix_free (rjX);
  gsl_matrix_free (rkX);
  gsl_matrix_free (E);
  gsl_matrix_free (temp);
}

/*
  calculate rotation diffusion tensor 

  Drr = kT [Krr - Ktr Ktt^-1 Krt]^-1	
  Dtt = kT [Ktt - Ktr Krr^-1 Krt]^-1	
  Dtr = -Ktt^-1 Ktr Drr
*/
void calc_diffusion_tensor (
  gsl_matrix *Ktt, gsl_matrix *Ktr, gsl_matrix *Krt, gsl_matrix *Krr,
  gsl_matrix *Dtt, gsl_matrix *Dtr, gsl_matrix *Drt, gsl_matrix *Drr)
{
  gsl_matrix *Ktt_LU  = gsl_matrix_alloc (3,3);
  gsl_matrix *Krr_LU  = gsl_matrix_alloc (3,3);
  gsl_matrix *Ktt_inv = gsl_matrix_alloc (3,3);
  gsl_matrix *Krr_inv = gsl_matrix_alloc (3,3);
  gsl_matrix *temp    = gsl_matrix_alloc (3,3);
  gsl_permutation *p  = gsl_permutation_alloc(3);
  int signum;

  gsl_matrix_memcpy(Ktt_LU,Ktt);
  gsl_matrix_memcpy(Krr_LU,Krr);
  gsl_permutation_init (p);
  gsl_linalg_LU_decomp (Ktt_LU,p,&signum);
  gsl_linalg_LU_invert (Ktt_LU,p,Ktt_inv);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ktr,Ktt_inv,0,temp);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,temp,Krt,1,Krr_LU);
  gsl_permutation_init (p);
  gsl_linalg_LU_decomp (Krr_LU,p,&signum);
  gsl_linalg_LU_invert (Krr_LU,p,Drr);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ktt_inv,Ktr,0,temp);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,temp,Drr,0,Dtr);
  gsl_matrix_memcpy(Ktt_LU,Ktt);
  gsl_matrix_memcpy(Krr_LU,Krr);
  gsl_permutation_init (p);
  gsl_linalg_LU_decomp (Krr_LU,p,&signum);
  gsl_linalg_LU_invert (Krr_LU,p,Krr_inv);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ktr,Krr_inv,0,temp);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,temp,Krt,1,Ktt_LU);
  gsl_permutation_init (p);
  gsl_linalg_LU_decomp (Ktt_LU,p,&signum);
  gsl_linalg_LU_invert (Ktt_LU,p,Dtt);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ktt_inv,Krt,0,temp);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,temp,Drr,0,Drt);
  //gsl_matrix_transpose_memcpy(Drt,Dtr);
  gsl_matrix_free (Ktt_LU);
  gsl_matrix_free (Krr_LU);
  gsl_matrix_free (Ktt_inv);
  gsl_matrix_free (Krr_inv);
  gsl_matrix_free (temp);
  gsl_permutation_free (p);
}

/*
  calculate center of diffusion
*/

void calc_center_of_diffusion (gsl_matrix *Drr,gsl_matrix *Dtr,gsl_vector *OR)
{
  gsl_matrix *CtDm     = gsl_matrix_alloc (3,3);
  gsl_matrix *CtDm_inv = gsl_matrix_alloc (3,3);
  gsl_vector *CtDv     = gsl_vector_alloc (3);
  gsl_permutation *p = gsl_permutation_alloc (3);
  int signum;

  gsl_matrix_set (CtDm,0,0, gsl_matrix_get(Drr,1,1)+gsl_matrix_get(Drr,2,2));
  gsl_matrix_set (CtDm,0,1,-gsl_matrix_get(Drr,0,1));
  gsl_matrix_set (CtDm,0,2,-gsl_matrix_get(Drr,0,2));
  gsl_matrix_set (CtDm,1,0,-gsl_matrix_get(Drr,0,1));
  gsl_matrix_set (CtDm,1,1, gsl_matrix_get(Drr,0,0)+gsl_matrix_get(Drr,2,2));
  gsl_matrix_set (CtDm,1,2,-gsl_matrix_get(Drr,1,2));
  gsl_matrix_set (CtDm,2,0,-gsl_matrix_get(Drr,0,2));
  gsl_matrix_set (CtDm,2,1,-gsl_matrix_get(Drr,1,2));
  gsl_matrix_set (CtDm,2,2, gsl_matrix_get(Drr,1,1)+gsl_matrix_get(Drr,0,0));

  gsl_permutation_init (p);
  gsl_linalg_LU_decomp (CtDm,p,&signum);
  gsl_linalg_LU_invert (CtDm,p,CtDm_inv);

  gsl_vector_set (CtDv,0, gsl_matrix_get(Dtr,1,2)-gsl_matrix_get(Dtr,2,1));
  gsl_vector_set (CtDv,1, gsl_matrix_get(Dtr,2,0)-gsl_matrix_get(Dtr,0,2));
  gsl_vector_set (CtDv,2, gsl_matrix_get(Dtr,0,1)-gsl_matrix_get(Dtr,1,0));
  gsl_blas_dgemv (CblasNoTrans,1,CtDm_inv,CtDv,0,OR);

  gsl_permutation_free (p);
  gsl_matrix_free (CtDm);
  gsl_matrix_free (CtDm_inv);
  gsl_vector_free (CtDv);
}

/* 	
  read MSMS .vert and .face files 
*/
size_t read_msms_surface (const char *prefix, gsl_matrix *f) 
{
  char line[MAXLENGTH],vert[MAXLENGTH],face[MAXLENGTH];
  size_t l,v1,v2,v3,vi,fi,nv,nf;
  double x,y,z;
  gsl_matrix *v;

  sprintf(vert,"%s.vert",prefix);
  sprintf(face,"%s.face",prefix);
  ifstream v_file (vert);
  ifstream f_file (face);
  if (f_file.fail() || v_file.fail()) {
    printf("error: cannot open %s and/or %s\n",vert,face);
    exit(-1);
    }

  if (f == NULL) {
    l = 0;
    while (! f_file.eof()) {
      f_file.getline(line,sizeof(line));
      l++;
      if (l == 3) {
	sscanf(line,"%d",&nf);
	break;
	}
      }
    f_file.close();
    return nf;
    }

  vi = l =  0;
  while (! v_file.eof()) {
    v_file.getline(line,sizeof(line));
    l++;
    if (l == 3) {
      sscanf(line,"%d",&nv);
      v = gsl_matrix_calloc (nv, 3);
      }
    if (l > 3 && strlen(line)) {
      sscanf(line,"%lf %lf %lf",&x,&y,&z);
      gsl_matrix_set (v, vi, 0, x);
      gsl_matrix_set (v, vi, 1, y);
      gsl_matrix_set (v, vi, 2, z);
      vi++;
      }
    }
  v_file.close();

  fi = l = 0;
  while (! f_file.eof()) {
    f_file.getline(line,sizeof(line));
    l++;
    if (l > 3 && strlen(line)) {
      sscanf(line,"%d %d %d",&v1,&v2,&v3);
      gsl_matrix_set (f, fi, 0, gsl_matrix_get(v, v1-1, 0));
      gsl_matrix_set (f, fi, 1, gsl_matrix_get(v, v1-1, 1));
      gsl_matrix_set (f, fi, 2, gsl_matrix_get(v, v1-1, 2));
      gsl_matrix_set (f, fi, 3, gsl_matrix_get(v, v2-1, 0));
      gsl_matrix_set (f, fi, 4, gsl_matrix_get(v, v2-1, 1));
      gsl_matrix_set (f, fi, 5, gsl_matrix_get(v, v2-1, 2));
      gsl_matrix_set (f, fi, 6, gsl_matrix_get(v, v3-1, 0));
      gsl_matrix_set (f, fi, 7, gsl_matrix_get(v, v3-1, 1));
      gsl_matrix_set (f, fi, 8, gsl_matrix_get(v, v3-1, 2));
      fi++;
      }
    }
  f_file.close();
  gsl_matrix_free (v);
  printf("MSMS_surface : %s %s (",vert,face);
  printf("Vertex: %d ",vi);
  printf("Face: %d)\n",fi);
  return nf;
}

/* 	
  read GTS surface 
*/

size_t read_gts_surface (const char *prefix, gsl_matrix *f) 
{
  char line[MAXLENGTH], gts[MAXLENGTH];
  size_t l,v1,v2,vi,ei,fi,nv,ne,nf;
  size_t i,j,k,imax,eidx[3],vidx[3];
  double x,y,z;
  gsl_matrix *v,*e;

  sprintf(gts,"%s.gts",prefix);
  ifstream gts_file (gts);
  if (gts_file.fail()) {
    printf("error: cannot open %s\n",gts);
    exit(-1);
    }

  if (f == NULL) {
    l = 0;
    while (! gts_file.eof()) {
      gts_file.getline(line,sizeof(line));
      l++;
      if (l == 1) {
	sscanf(line,"%d %d %d",&nv,&ne,&nf);
        break;
        }
      }
    gts_file.close();
    return nf;
    }

  l = 0;
  while (! gts_file.eof()) {
    gts_file.getline(line,sizeof(line));
    l++;
    if (l == 1) {
      sscanf(line,"%d %d %d",&nv,&ne,&nf);
      v = gsl_matrix_alloc (nv, 3);
      e = gsl_matrix_alloc (ne, 2);
      vi = ei = fi = 0;
      }
    /* vertices */
    if (l > 1 && l <= (nv+1) && strlen(line)) {
      sscanf(line,"%lf %lf %lf",&x,&y,&z);
      gsl_matrix_set (v, vi, 0, x);
      gsl_matrix_set (v, vi, 1, y);
      gsl_matrix_set (v, vi, 2, z);
      vi++;
      }
    /* edges */
    if (l > (nv+1) && l <= (nv+ne+1) && strlen(line)) {
      sscanf(line,"%d %d",&v1,&v2);
      gsl_matrix_set (e, ei, 0, v1-1);
      gsl_matrix_set (e, ei, 1, v2-1);
      ei++;
      }
    /* faces */
    if (l > (nv+ne+1) && strlen(line)) {
      sscanf(line,"%d %d %d",&eidx[0],&eidx[1],&eidx[2]);
      for(imax=0,k=0;k<3;k++) {
	for (j=0;j<2;j++) {
	  for (i=0;i<imax;i++)
	    if (vidx[i] == gsl_matrix_get(e,eidx[k]-1,j)) break;
	  if (i == imax) 
	    vidx[imax++] = gsl_matrix_get(e,eidx[k]-1,j);
	  }
	}
      gsl_matrix_set (f, fi, 0, gsl_matrix_get(v, vidx[0], 0));
      gsl_matrix_set (f, fi, 1, gsl_matrix_get(v, vidx[0], 1));
      gsl_matrix_set (f, fi, 2, gsl_matrix_get(v, vidx[0], 2));
      gsl_matrix_set (f, fi, 3, gsl_matrix_get(v, vidx[1], 0));
      gsl_matrix_set (f, fi, 4, gsl_matrix_get(v, vidx[1], 1));
      gsl_matrix_set (f, fi, 5, gsl_matrix_get(v, vidx[1], 2));
      gsl_matrix_set (f, fi, 6, gsl_matrix_get(v, vidx[2], 0));
      gsl_matrix_set (f, fi, 7, gsl_matrix_get(v, vidx[2], 1));
      gsl_matrix_set (f, fi, 8, gsl_matrix_get(v, vidx[2], 2));
      fi++;
      }
    }
  gts_file.close();
  gsl_matrix_free (v);
  gsl_matrix_free (e);
  printf("GTS_surface  : %s (",gts);
  printf("Vertex: %d ",vi);
  printf("Edge: %d ",ei);
  printf("Face: %d)\n",fi);
  return nf;
}

/* 	
  read GeomView surface  (OFF format)
*/

size_t read_geom_surface (const char *prefix, gsl_matrix *f) 
{
  char line[MAXLENGTH],dummy[MAXLENGTH],off[MAXLENGTH];
  size_t l,v1,v2,v3,vi,fi,nv,ne,nf,dim;
  double x,y,z;
  gsl_matrix *v;

  sprintf(off,"%s.off",prefix);
  ifstream geom_file (off);
  if (geom_file.fail()) {
    printf("error: cannot open %s\n",off);
    exit (-1);
    }

  if (f == NULL) {
    l = 0;
    while (! geom_file.eof()) {
      geom_file.getline(line,sizeof(line));
      l++;
      if (l == 1) {
	sscanf(line,"%s %d %d %d",dummy,&nv,&nf,&ne);
	break;
	}
      }
    geom_file.close();
    return nf;
    }

  l = 0;
  while (! geom_file.eof()) {
    geom_file.getline(line,sizeof(line));
    l++;
    if (l == 1) {
      sscanf(line,"%s %d %d %d",dummy,&nv,&nf,&ne);
      v = gsl_matrix_alloc (nv, 3);
      vi = fi = 0;
      }
    if (l > 1 && l <= (nv+1) && strlen(line)) {
      sscanf(line,"%lf %lf %lf",&x,&y,&z);
      gsl_matrix_set (v, vi, 0, x);
      gsl_matrix_set (v, vi, 1, y);
      gsl_matrix_set (v, vi, 2, z);
      vi++;
      }
    if (l > (nv+1) && strlen(line)) {
      sscanf(line,"%d %d %d %d",&dim,&v1,&v2,&v3);
      gsl_matrix_set (f, fi, 0, gsl_matrix_get(v, v1, 0));
      gsl_matrix_set (f, fi, 1, gsl_matrix_get(v, v1, 1));
      gsl_matrix_set (f, fi, 2, gsl_matrix_get(v, v1, 2));
      gsl_matrix_set (f, fi, 3, gsl_matrix_get(v, v2, 0));
      gsl_matrix_set (f, fi, 4, gsl_matrix_get(v, v2, 1));
      gsl_matrix_set (f, fi, 5, gsl_matrix_get(v, v2, 2));
      gsl_matrix_set (f, fi, 6, gsl_matrix_get(v, v3, 0));
      gsl_matrix_set (f, fi, 7, gsl_matrix_get(v, v3, 1));
      gsl_matrix_set (f, fi, 8, gsl_matrix_get(v, v3, 2));
      fi++;
      }
    }
  geom_file.close();
  gsl_matrix_free (v);
  printf("GEOM_surface : %s (",off);
  printf("Vertex: %d ",vi);
  printf("Face: %d)\n",fi);
  return nf;
}

/* 	
  write GeomView surface with face color (OFF format)
*/

void write_geom_surface (const char *prefix, gsl_matrix *f, gsl_vector *vf,
  int color_scheme) 
{
  char off[MAXLENGTH];
  size_t N,fi;
  double color_r,color_g,color_b,color_a;

  sprintf(off,"%s.%d-%d.off",prefix,(int)gboundary,(int)epsilon);
  ofstream geom_file (off);
  if (geom_file.fail()) {
    printf("error: cannot open %s\n",off);
    exit (-1);
    }

  if (f != NULL) {
    N = f->size1;
    geom_file << "OFF\n";
    geom_file << 3*N << " " << N << " " << 3*N/2 << "\n";
    fi = 0;
    /* vertices */
    while (fi < N) {
      geom_file << gsl_matrix_get(f,fi,0) << " ";
      geom_file << gsl_matrix_get(f,fi,1) << " ";
      geom_file << gsl_matrix_get(f,fi,2) << "\n";
      geom_file << gsl_matrix_get(f,fi,3) << " ";
      geom_file << gsl_matrix_get(f,fi,4) << " ";
      geom_file << gsl_matrix_get(f,fi,5) << "\n";
      geom_file << gsl_matrix_get(f,fi,6) << " ";
      geom_file << gsl_matrix_get(f,fi,7) << " ";
      geom_file << gsl_matrix_get(f,fi,8) << "\n";
      fi++;
      }
    /* faces */
    fi = 0;
    while (fi < N) {
      /* vertices indices */
      geom_file << "3 "<< 3*fi << " " << 3*fi+1 << " " << 3*fi+2 << " ";
      if (color_scheme == 0) {
	/* 
	  RGB colormap according to velocity field correlation factor
	  vfcf=1 blue, 
	  vfcf=0 red 
	*/
	color_r = 1.0-gsl_vector_get(vf,fi);
	color_g = 0.0;
	color_b = gsl_vector_get(vf,fi);
	color_a = 0.75;
	}
      if (color_scheme == 1) {
	/* 
	  RGB colormap according to N|C terminal area 
	  N-terminus: blue
	  C-terminus: red
	  Other: white
	*/
	color_r = 1.0;
	color_g = 1.0-gsl_vector_get(vf,fi);
	color_b = 1.0-gsl_vector_get(vf,fi);
	color_a = 0.75;
	}
      geom_file << color_r << " " << color_g << " " << color_b << " ";
      geom_file << color_a << "\n"; /* opacity or transparency */
      fi++;
      }
    }
  geom_file.close();
}

size_t read_pdb_coor (const char *prefix, vector <Atom> & coor)
{
  char filename[MAXLENGTH];
  char line[MAXLENGTH];
  string s,atmName,resNum,resName,atmCoor,chainId;
  unsigned int c;
  Atom a;

  sprintf(filename,"%s.pdb",prefix);
  ifstream pdb_file (filename);

  if (pdb_file.fail()) {
    printf("error: cannot open %s\n",filename);
    exit (-1);
    }
  
  coor.resize(0);

  while (! pdb_file.eof()) {
    pdb_file.getline(line, sizeof(line));
    if (strlen(line) > 0) {
      s = line;
      if(s.find("ATOM") != string::npos || s.find("HETATM") != string::npos) {
	// PDB format version 2.3
	atmName = s.substr(12,4);
	resName = s.substr(17,3);
	chainId = s.substr(21,1);
	resNum  = s.substr(22,4);
	atmCoor = s.substr(30,24);
	while ((c=atmName.find(" ",0)) != string::npos)
	  atmName.replace(c,1,"");// trim empty spaces
	a.name = atmName;
	a.resn = resName;
	sscanf(chainId.c_str(),"%c",&a.c);
	sscanf(resNum.c_str(),"%d",&a.r);
        sscanf(atmCoor.c_str(),"%lf %lf %lf",&a.x,&a.y,&a.z);
	coor.push_back(a);
	} // ATOM or HETATM
      } // if line is not empty
    } // until end of file

  return coor.size();

}

void time_elapsed (time_t &a,time_t &b) 
{
  int hh=0,mm=0,ss=0,elapsed=b-a;
  hh = (int)floor((double)(elapsed/3600));
  mm = (int)floor((double)((elapsed/60)-(hh*60)));
  ss = (elapsed - (hh*3600)-(mm*60));
  printf("time %02d:%02d:%02d\n",hh,mm,ss);
}

/* 
  GSL functions for numerical integration 
  Numerical Integration by QUADPACK through GSL	
    |RESULT − I| ≤ max(epsabs, epsrel |I|) 
    The QAG algorithm is a simple adaptive integration procedure.
    The integration region is divided into subintervals, 
    and on each iteration the subinterval with the largest estimated error
    is bisected. This reduces the overall error rapidly, as the subintervals
    become concentrated around local difficulties in the integrand.
    These subintervals are managed by a gsl_integration_workspace struct, 
    which handles the memory for the subinterval ranges, results and 
    error estimates.
*/

double dyt0 (double y, void * params) 
{
  scalar_product s = *(scalar_product *) params;
  double n = s.ad2+y*s.ac+s.a*sqrt(y*y*s.sqc+2*y*s.cd2+s.sqd2);
  double d = s.ad1+y*s.ab+s.a*sqrt(y*y*s.sqb+2*y*s.bd1+s.sqd1);
  return (1/s.a)*log(n/d);
}
double dyg0 (double y, void * params) 
{
  scalar_product s = *(scalar_product *) params;
  double n1 = (s.ad1+y*s.ab)/sqrt(s.sqd1+2*y*s.bd1+y*y*s.sqb);
  double n2 = (s.ad2+y*s.ac)/sqrt(s.sqd2+2*y*s.cd2+y*y*s.sqc);
  double d1 = s.ad1*s.ad1-s.sqd1*s.sqa-2*y*s.bd1*s.sqa;
  double d2 = 2*y*s.ad1*s.ab+y*y*s.ab*s.ab-y*y*s.sqa*s.sqb;
  return (n1-n2)/(d1+d2);
}
double dyyg0 (double y, void * params) 
{
  scalar_product s = *(scalar_product *) params;
  double n1 = (s.ad1+y*s.ab)/sqrt(s.sqd1+2*y*s.bd1+y*y*s.sqb);
  double n2 = (s.ad2+y*s.ac)/sqrt(s.sqd2+2*y*s.cd2+y*y*s.sqc);
  double d1 = s.ad1*s.ad1-s.sqd1*s.sqa-2*y*s.bd1*s.sqa;
  double d2 = 2*y*s.ad1*s.ab+y*y*s.ab*s.ab-y*y*s.sqa*s.sqb;
  return y*(n1-n2)/(d1+d2);
}
double dyyyg0 (double y, void * params) 
{
  scalar_product s = *(scalar_product *) params;
  double n1 = (s.ad1+y*s.ab)/sqrt(s.sqd1+2*y*s.bd1+y*y*s.sqb);
  double n2 = (s.ad2+y*s.ac)/sqrt(s.sqd2+2*y*s.cd2+y*y*s.sqc);
  double d1 = s.ad1*s.ad1-s.sqd1*s.sqa-2*y*s.bd1*s.sqa;
  double d2 = 2*y*s.ad1*s.ab+y*y*s.ab*s.ab-y*y*s.sqa*s.sqb;
  return y*y*(n1-n2)/(d1+d2);
}
double dyg1 (double y, void * params) 
{
  scalar_product s = *(scalar_product *) params;
  double n1 = sqrt(s.sqd1+2*y*s.bd1+y*y*s.sqb);
  double n2 = s.d1d2-y*(s.ad1-2*s.bd1-s.ab)+y*y*s.bc;
  double n3 = sqrt(s.sqd2+2*y*s.cd2+y*y*s.sqc);
  double d1 = s.ad1*s.ad1-s.sqd1*s.sqa-2*y*s.bd1*s.sqa;
  double d2 = 2*y*s.ad1*s.ab+y*y*s.ab*s.ab-y*y*s.sqa*s.sqb;
  return (-n1+n2/n3)/(d1+d2);
}
double dyyg1 (double y, void * params) 
{
  scalar_product s = *(scalar_product *) params;
  double n1 = sqrt(s.sqd1+2*y*s.bd1+y*y*s.sqb);
  double n2 = s.d1d2-y*(s.ad1-2*s.bd1-s.ab)+y*y*s.bc;
  double n3 = sqrt(s.sqd2+2*y*s.cd2+y*y*s.sqc);
  double d1 = s.ad1*s.ad1-s.sqd1*s.sqa-2*y*s.bd1*s.sqa;
  double d2 = 2*y*s.ad1*s.ab+y*y*s.ab*s.ab-y*y*s.sqa*s.sqb;
  return y*(-n1+n2/n3)/(d1+d2);
}
double dxxxg2 (double x, void * params) 
{
  scalar_product s = *(scalar_product *) params;
  double n1 = (s.bd1+x*s.ab)/sqrt(s.sqd1+2*x*s.ad1+x*x*s.sqa);
  double n2 = (s.bd3-x*s.bc)/sqrt(s.sqd3-2*x*s.cd3+x*x*s.sqc);
  double d1 = s.bd1*s.bd1-s.sqb*s.sqd1-2*x*s.ad1*s.sqb;
  double d2 = 2*x*s.bd1*s.ab+x*x*s.ab*s.ab-x*x*s.sqa*s.sqb;
  return x*x*(n1-n2)/(d1+d2);
}

void calc_surf_geometry (
    gsl_matrix *Vert, 
    gsl_matrix *Centroid,
    gsl_matrix *Side, 
    gsl_matrix *Length,
    gsl_vector *Area)
{
  size_t k;
  double a,b,c,s;
  gsl_vector_view v1,v2,v3,va,vb,vc,vk,sk,ck;

  for (k=0;k<Vert->size1;k++) {

    /* vertices */ 
    vk = gsl_matrix_row (Vert,k);	
    v1 = gsl_vector_subvector (&vk.vector,0,3);
    v2 = gsl_vector_subvector (&vk.vector,3,3);
    v3 = gsl_vector_subvector (&vk.vector,6,3);

    /* sides */
    if (Side != NULL && Length != NULL && Area != NULL) {
      sk = gsl_matrix_row (Side,k);
      va = gsl_vector_subvector (&sk.vector,0,3);
      vb = gsl_vector_subvector (&sk.vector,3,3);
      vc = gsl_vector_subvector (&sk.vector,6,3);
      /* k, vector va = v2 - v1, a */
      gsl_vector_memcpy(&va.vector,&v2.vector);
      gsl_vector_sub(&va.vector,&v1.vector);
      a = gsl_blas_dnrm2(&va.vector);
      gsl_matrix_set(Length,k,0,a);
      /* k, vector vb = v3 - v1, b */
      gsl_vector_memcpy(&vb.vector,&v3.vector);
      gsl_vector_sub(&vb.vector,&v1.vector);
      b = gsl_blas_dnrm2(&vb.vector);
      gsl_matrix_set(Length,k,1,b);
      /* k, vector vc = v3 - v2, c */
      gsl_vector_memcpy(&vc.vector,&v3.vector);
      gsl_vector_sub(&vc.vector,&v2.vector);
      c = gsl_blas_dnrm2(&vc.vector);
      gsl_matrix_set(Length,k,2,c);
      /* surface area using Heron's formula */
      s = 0.5*(a+b+c); /* semiperimeter */
      gsl_vector_set (Area,k, sqrt(s*(s-a)*(s-b)*(s-c)));
      }
  
    /* k, centroid of triangle */
    if (Centroid != NULL) {
      ck = gsl_matrix_row (Centroid,k);
      gsl_vector_set_zero (&ck.vector);
      gsl_vector_add (&ck.vector,&v3.vector);
      gsl_vector_add (&ck.vector,&v2.vector);
      gsl_vector_add (&ck.vector,&v1.vector);
      gsl_vector_scale (&ck.vector,(1.0/3.0));
      }

    }// k
}

void calc_velocity_factor (gsl_matrix *C,Attrib & dSa,gsl_vector *vf)
{
  extern double epsilon, gboundary;

  gsl_histogram *g = gsl_histogram_alloc (29);
  gsl_histogram *h = gsl_histogram_alloc (10);
  double range_g [30] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,
			3.0,3.3,3.6,3.9,4.2,4.5,4.8,5.1,5.4,5.7,
			6.0,6.3,6.9,7.2,7.5,7.8,8.1,8.4,8.7,1000.};
  double range_h[11] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1};
  gsl_histogram_set_ranges (g, range_g, 30);
  gsl_histogram_set_ranges (h, range_h, 11);

  double d,dmin,lower,upper;
  size_t k;

  for (k=0;k<C->size1;k++) {
    dmin = dSa.d[k];
    d = _max (dmin - gboundary, 0.);
    gsl_vector_set (vf,k,exp(-d/epsilon)); // velocity_factor = exp(-d/epsilon)
    if (verb) {
      gsl_histogram_increment (g,dmin);
      gsl_histogram_increment (h,exp(-d/epsilon));
      } // verb
    } // k

  /* histogram */
  if (verb) {
    printf("Offset distance \t");
    printf("Velocity_Factor \n");
    for (k=0;k<g->n;k++) {
      if (k == g->n-1) {
	gsl_histogram_get_range (g,k,&lower,&upper);
	printf("%3.1f <     %5g\t\t",lower,gsl_histogram_get (g,k));
	}
      else {
	gsl_histogram_get_range (g,k,&lower,&upper);
	printf("%3.1f - %3.1f %5g\t\t",lower,upper,gsl_histogram_get (g,k));
	}
      if (k < h->n) {
	if (k == 0) {
	  gsl_histogram_get_range (h,h->n-k-1,&lower,&upper);
	  printf("%3.1f - 1.0 %5g\n",lower,gsl_histogram_get (h,h->n-k-1));
	  }
	else {
	  gsl_histogram_get_range (h,h->n-k-1,&lower,&upper);
	  printf("%3.1f - %3.1f %5g\n",lower,upper,
	    gsl_histogram_get (h,h->n-k-1));
	  }
	}
      else printf("\n");
      }
    }// verb

  gsl_histogram_free (g);
  gsl_histogram_free (h);
}

void calc_dimension (gsl_matrix *Vert) 
{
  gsl_vector_view v1,v2,v3;
  double dim_x,dim_y,dim_z;

  v1 = gsl_matrix_column (Vert,0); /* v1_x */
  v2 = gsl_matrix_column (Vert,3); /* v2_x */
  v3 = gsl_matrix_column (Vert,6); /* v3_x */
  dim_x = 
    _max (
      _max (gsl_vector_max(&v1.vector),gsl_vector_max(&v2.vector)),
      gsl_vector_max(&v3.vector))-
    _min (
      _min (gsl_vector_min(&v1.vector),gsl_vector_min(&v2.vector)),
      gsl_vector_min(&v3.vector));
  v1 = gsl_matrix_column (Vert,1); /* v1_y */
  v2 = gsl_matrix_column (Vert,4); /* v2_y */
  v3 = gsl_matrix_column (Vert,7); /* v3_y */
  dim_y = 
    _max (
      _max (gsl_vector_max(&v1.vector),gsl_vector_max(&v2.vector)),
      gsl_vector_max(&v3.vector))-
    _min (
      _min (gsl_vector_min(&v1.vector),gsl_vector_min(&v2.vector)),
      gsl_vector_min(&v3.vector));
  v1 = gsl_matrix_column (Vert,2); /* v1_z */
  v2 = gsl_matrix_column (Vert,5); /* v2_z */
  v3 = gsl_matrix_column (Vert,8); /* v3_z */
  dim_z = 
    _max (
      _max (gsl_vector_max(&v1.vector),gsl_vector_max(&v2.vector)),
      gsl_vector_max(&v3.vector))-
    _min (
      _min (gsl_vector_min(&v1.vector),gsl_vector_min(&v2.vector)),
      gsl_vector_min(&v3.vector));

  printf("Dimension    : X %9.4f  Y %9.4f Z %9.4f\n",dim_x,dim_y,dim_z);
}

void integ_err (size_t k,size_t j,scalar_product &D) 
{
  fprintf(stderr,"integration error at k=%d j=%d\n",k,j);
  fprintf(stderr,"ad1:%g ad2:%g ac:%g ab:%g bc:%g\n",
    D.ad1,D.ad2,D.ac,D.ab,D.bc);
  fprintf(stderr,"bd1:%g bd3:%g cd2:%g cd3:%g d1d2:%g\n",
    D.bd1,D.bd3,D.cd2,D.cd3,D.d1d2);
  fprintf(stderr,"a:%g a*a:%g b*b:%g c*c:%g\n",
    D.a,D.sqa,D.sqb,D.sqc);
  fprintf(stderr,"d1*d1:%g d2*d2:%g d3*d3:%g\n",
    D.sqd1,D.sqd2,D.sqd3);
  exit(-1);
}

void calc_attrib (gsl_matrix *C, gsl_matrix *RC,vector <Atom> & coor,
  Attrib &dSa)
{
  extern int algorithm;
  
  dSa.atom.resize(0); 
  dSa.resid.resize(0);
  dSa.chain.resize(0);
  dSa.anchor_idx.resize(0);
  dSa.anchor_element.resize(0);
  dSa.delta_r.resize(0);
  dSa.d.resize(0);

  gsl_vector *dist = gsl_vector_alloc (3);
  gsl_vector_view ck,cr;
  size_t k,r;
  double d, dmin;

  /* 
  algorithms for d calculation:

  for a given soft surface element,

  [1] search for a closest surface element within the rigid surface elements
      use a distance to the centroid (default)
  [2] search for a rigid/soft surfaces junction element(s)
      use a distance to the junction element(s)
  [3] search for a rigid/soft surfaces junction residue
      use a distance to the junction residue along the backbone covalent bonds
  */

  //
  // ALGORITHM 1 (default)
  //
  if (algorithm == 1) {
    for (k=0;k<C->size1;k++) {
      ck = gsl_matrix_row (C,k);
      dmin = 1.0e+6;
      for (r=0;r<RC->size1;r++) {
	cr = gsl_matrix_row (RC,r);
	gsl_vector_memcpy (dist,&ck.vector);
	gsl_vector_sub (dist,&cr.vector); 
	d  = gsl_blas_dnrm2(dist); /* d = |C_k - RC_r| */
	if (dmin > d) dmin = d;
	} // r
      dSa.d.push_back (dmin);
      } // k
    return;
    } // algorithm 1 or method (a)

  extern vector < char > Rchain;
  extern vector < int > RNter,RCter;

  extern double assoc_threshold;
  double assoc_threshold_inc;

  gsl_vector_view ci;
  gsl_vector *cp = gsl_vector_alloc (3);
  vector <int> atomlist, residlist, residfreq;
  vector <int> anchor_n_list, anchor_c_list, anchor_NULL;
  vector <char> chainlist;
  char chain;
  int resid, maxfreq, anchor_n, anchor_c;
  size_t i,c,p,neighbor;

  // COMMON for ALGORITHMS 2 and 3
  // connect an element to a specific residue
  for (k=0;k<C->size1;k++) {
    ck = gsl_matrix_row (C,k);
    atomlist.resize(0);
    residlist.resize(0);
    residfreq.resize(0);
    assoc_threshold_inc = 0.;
    while (residlist.size() == 0) {
      for (p=0;p<coor.size();p++) {
	gsl_vector_memcpy (dist,&ck.vector);
	gsl_vector_set (cp, 0, coor[p].x);
	gsl_vector_set (cp, 1, coor[p].y);
	gsl_vector_set (cp, 2, coor[p].z);
	gsl_vector_sub (dist, cp);

	d = gsl_blas_dnrm2 (dist);

	if (d < assoc_threshold + assoc_threshold_inc) {
	  atomlist.push_back (p);
	  /* make a non-redundant resid list */
	  for (i=0;i<residlist.size();i++) {
	    if (residlist[i] == coor[p].r) break;
	    }
	  if (i == residlist.size()) {
	    chainlist.push_back(coor[p].c);
	    residlist.push_back(coor[p].r);
	    residfreq.push_back(1);
	    }
	  else {
	    residfreq[i] += 1;
	    }
	  } // if d < assoc_threshold

	} // for p


      if (residlist.size() > 0) {
	/* choose most relevant residue */
	chain = 'A';
	resid = -1;
	maxfreq = 0;
	for (i=0;i<residlist.size();i++) {
	  if (residfreq[i] > maxfreq) {
	    maxfreq = residfreq[i];
	    chain = chainlist[i];
	    resid = residlist[i];
	    }
	  }
	dSa.atom.push_back (atomlist);
	dSa.chain.push_back (chain);
	dSa.resid.push_back (resid);
	}
      else {
	assoc_threshold_inc += 2.;
	}
      } // repeat unless relevant residue is found
    } // k

  /* anchor point and residue distance */
  for (c=0;c<Rchain.size();c++) {

    anchor_n_list.resize(0);
    anchor_c_list.resize(0);
    anchor_NULL.resize(0);

    if (algorithm == 2) {
      // find N-terminus element(s) within the reference surface elements
      assoc_threshold_inc = 0.;
      while (anchor_n_list.size() == 0)  { // ensure at leat one element
	for (p=0;p<coor.size();p++) {
	  if (coor[p].c == Rchain[c] && coor[p].r == RNter[c]) {
	    gsl_vector_set (cp, 0, coor[p].x);
	    gsl_vector_set (cp, 1, coor[p].y);
	    gsl_vector_set (cp, 2, coor[p].z);
	    for (r=0;r<RC->size1;r++) {
	      cr = gsl_matrix_row (RC, r);
	      gsl_vector_memcpy (dist, &cr.vector);
	      gsl_vector_sub (dist, cp);
	      d = gsl_blas_dnrm2 (dist);
	      if (d < assoc_threshold + assoc_threshold_inc)
		anchor_n_list.push_back (r);
	      }
	    } // = RNter
	  } // p, coor.size()
	assoc_threshold_inc += 2.;
	} // while
      // find C-terminus element(s) within the reference surface elements
      assoc_threshold_inc = 0.;
      while (anchor_c_list.size() == 0)  { // ensure at least one element
	for (p=0;p<coor.size();p++) {
	  if (coor[p].c == Rchain[c] && coor[p].r == RCter[c]) {
	    gsl_vector_set (cp, 0, coor[p].x);
	    gsl_vector_set (cp, 1, coor[p].y);
	    gsl_vector_set (cp, 2, coor[p].z);
	    for (r=0;r<RC->size1;r++) {
	      cr = gsl_matrix_row (RC, r);
	      gsl_vector_memcpy (dist, &cr.vector);
	      gsl_vector_sub (dist, cp);
	      d = gsl_blas_dnrm2 (dist);
	      if (d < assoc_threshold + assoc_threshold_inc)
		anchor_c_list.push_back (r);
	      }
	    } // = RCter
	  } // p, coor.size()
	assoc_threshold_inc += 2.;
	} // while
      // store elements in anchor_element
      for (k=0;k<C->size1;k++) {
	if (dSa.chain[k] == Rchain[c]) {
	  if (dSa.resid[k] <= RNter[c])
	    dSa.anchor_element.push_back (anchor_n_list);
	  else if (dSa.resid[k] >= RCter[c])
	    dSa.anchor_element.push_back (anchor_c_list);
	  else
	    dSa.anchor_element.push_back (anchor_NULL);
	  }
	} // k
      }// algorithm-2 or method (c)

    //
    // ALGORITHM 3 or method (d)
    //
    if (algorithm == 3) {
      for (k=0;k<C->size1;k++) {
	if (dSa.chain[k] == Rchain[c]) {
	  if (dSa.resid[k] <= RNter[c])
	    dSa.delta_r.push_back (_abs (RNter[c] - dSa.resid[k]));
	  else if (dSa.resid[k] >= RCter[c])
	    dSa.delta_r.push_back (_abs (dSa.resid[k] - RCter[c]));
	  else
	    dSa.delta_r.push_back (-1); // not applicable
	  }
	} // k
      } // algorithm-3 or method (d)

    } // c, Rchain

  gsl_vector_free (cp);

  //
  // ALGORITHM 2 
  // 
  if (algorithm == 2) {
    gsl_vector *anchor = gsl_vector_alloc (3);
    for (k=0;k<C->size1;k++) {
      gsl_vector_set_zero (anchor);
      ck = gsl_matrix_row (C,k);
      if (dSa.anchor_element[k].size() > 0) {
	// average position of N- and C- terminus associated elements
	for (r=0;r<dSa.anchor_element[k].size();r++) {
	  cr = gsl_matrix_row (RC,dSa.anchor_element[k][r]);
	  gsl_blas_daxpy (1,&cr.vector,anchor);
	  } // r
	gsl_vector_scale (anchor,1./r);
	gsl_vector_memcpy (dist,&ck.vector);
        gsl_vector_sub (dist, anchor);
        dmin = gsl_blas_dnrm2 (dist); /* dmin = |C_k - ON| */
        }
      else {
	dmin = 1.0e+6;
	for (r=0;r<RC->size1;r++) {
	  cr = gsl_matrix_row (RC,r);
	  gsl_vector_memcpy (dist,&ck.vector);
	  gsl_vector_sub (dist,&cr.vector); 
	  d  = gsl_blas_dnrm2(dist); /* d = |C_k - RC_r| */
	  if (dmin > d) dmin = d;
	  } // r
        }
      dSa.d.push_back (dmin);
      } // k
    gsl_vector_free (anchor);
    } // algorithm 2

  if (algorithm == 3) {
    for (k=0;k<C->size1;k++)
      dSa.d.push_back (_max (dSa.delta_r[k], 0.));
    } // algorithm 3

}

/*
  Block Inversion
  
    M^-1  = | A B | ^-1 =  | O P |
	    | C D |        | Q R |

    O = A^-1 + (A^-1) B ((D-C(A^-1)B)^-1) C A^-1
    P = -(A^-1) B ((D-C(A^-1)B)^-1)
    Q = -(D-C(A^-1)B)^-1) C A^-1
    R = (D-C(A^-1)B)^-1)
    
    note: D-C(A^-1)B is called Schur complement of A
    
    (row_position,column_position) : row_size x column_size

    Matrix M : (0,0) (m+n) x (m+n)

    Matrix A : (0,0) m x m
    Matrix B : (0,m) m x n
    Matrix C : (m,0) n x m
    Matrix D : (m,m) n x n

    note: m and n are arbitrary size
*/

void block_inv (gsl_matrix *G, gsl_matrix *G_inv)
{
  size_t m = (size_t)(G->size1*0.5); 
  size_t n = G->size1 - m;
  gsl_matrix *block_A = gsl_matrix_alloc (m,m); // block A 
  gsl_matrix *Schur_A = gsl_matrix_alloc (n,n); // Schur complement of A
  gsl_matrix_view GA,GB,GC,GD,GiA,GiB,GiC,GiD;

  GA  = gsl_matrix_submatrix (G,0,0,m,m);
  GB  = gsl_matrix_submatrix (G,0,m,m,n);
  GC  = gsl_matrix_submatrix (G,m,0,n,m);
  GD  = gsl_matrix_submatrix (G,m,m,n,n);
  GiA = gsl_matrix_submatrix (G_inv,0,0,m,m);
  GiB = gsl_matrix_submatrix (G_inv,0,m,m,n);
  GiC = gsl_matrix_submatrix (G_inv,m,0,n,m);
  GiD = gsl_matrix_submatrix (G_inv,m,m,n,n);

  // block A and Schur complement of A
  gsl_matrix_memcpy (block_A, &GA.matrix);
  gsl_matrix_memcpy (Schur_A, &GD.matrix);

  // inversion of block A
  // = block A of inverse matrix (temporary)
  block_inv_2 (block_A, &GiA.matrix);

  // Ai_B (mxn) = A_inv (mxm) B (mxn)
  // C_Ai (nxm) = C (nxm) A_inv (mxm)
  gsl_matrix *Ai_B = gsl_matrix_alloc (m,n);
  gsl_matrix *C_Ai = gsl_matrix_alloc (n,m);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GiA.matrix,&GB.matrix,0,Ai_B);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GC.matrix,&GiA.matrix,0,C_Ai);

  // Shcur complement of block A (nxn) = C (nxm) A_inv (mxm) B (mxn)
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,C_Ai,&GB.matrix,1,Schur_A);

  // inversion of Shcur complement of block A (nxn) 
  // = block D of inverse matrix 
  block_inv_2 (Schur_A, &GiD.matrix);

  // block C of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,&GiD.matrix,C_Ai,0,&GiC.matrix);
  // block B of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiD.matrix,0,&GiB.matrix);
  // block A of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiC.matrix,1,&GiA.matrix);

  gsl_matrix_free (block_A);
  gsl_matrix_free (Schur_A);
  gsl_matrix_free (C_Ai);
  gsl_matrix_free (Ai_B);
}

void block_inv_2 (gsl_matrix *G, gsl_matrix *G_inv)
{
  size_t m = (size_t)(G->size1*0.5); 
  size_t n = G->size1 - m;
  gsl_matrix *block_A = gsl_matrix_alloc (m,m); // block A 
  gsl_matrix *Schur_A = gsl_matrix_alloc (n,n); // Schur complement of A
  gsl_matrix_view GA,GB,GC,GD,GiA,GiB,GiC,GiD;

  GA  = gsl_matrix_submatrix (G,0,0,m,m);
  GB  = gsl_matrix_submatrix (G,0,m,m,n);
  GC  = gsl_matrix_submatrix (G,m,0,n,m);
  GD  = gsl_matrix_submatrix (G,m,m,n,n);
  GiA = gsl_matrix_submatrix (G_inv,0,0,m,m);
  GiB = gsl_matrix_submatrix (G_inv,0,m,m,n);
  GiC = gsl_matrix_submatrix (G_inv,m,0,n,m);
  GiD = gsl_matrix_submatrix (G_inv,m,m,n,n);

  // block A and Schur complement of A
  gsl_matrix_memcpy (block_A, &GA.matrix);
  gsl_matrix_memcpy (Schur_A, &GD.matrix);

  // inversion of block A
  // = block A of inverse matrix (temporary)
  block_inv_3 (block_A, &GiA.matrix);

  // Ai_B (mxn) = A_inv (mxm) B (mxn)
  // C_Ai (nxm) = C (nxm) A_inv (mxm)
  gsl_matrix *Ai_B = gsl_matrix_alloc (m,n);
  gsl_matrix *C_Ai = gsl_matrix_alloc (n,m);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GiA.matrix,&GB.matrix,0,Ai_B);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GC.matrix,&GiA.matrix,0,C_Ai);

  // Shcur complement of block A (nxn) = C (nxm) A_inv (mxm) B (mxn)
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,C_Ai,&GB.matrix,1,Schur_A);

  // inversion of Shcur complement of block A (nxn) 
  // = block D of inverse matrix 
  block_inv_3 (Schur_A, &GiD.matrix);

  // block C of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,&GiD.matrix,C_Ai,0,&GiC.matrix);
  // block B of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiD.matrix,0,&GiB.matrix);
  // block A of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiC.matrix,1,&GiA.matrix);

  gsl_matrix_free (block_A);
  gsl_matrix_free (Schur_A);
  gsl_matrix_free (C_Ai);
  gsl_matrix_free (Ai_B);
}

void block_inv_3 (gsl_matrix *G, gsl_matrix *G_inv)
{
  size_t m = (size_t)(G->size1*0.5); 
  size_t n = G->size1 - m;
  gsl_matrix *block_A = gsl_matrix_alloc (m,m); // block A 
  gsl_matrix *Schur_A = gsl_matrix_alloc (n,n); // Schur complement of A
  gsl_matrix_view GA,GB,GC,GD,GiA,GiB,GiC,GiD;

  GA  = gsl_matrix_submatrix (G,0,0,m,m);
  GB  = gsl_matrix_submatrix (G,0,m,m,n);
  GC  = gsl_matrix_submatrix (G,m,0,n,m);
  GD  = gsl_matrix_submatrix (G,m,m,n,n);
  GiA = gsl_matrix_submatrix (G_inv,0,0,m,m);
  GiB = gsl_matrix_submatrix (G_inv,0,m,m,n);
  GiC = gsl_matrix_submatrix (G_inv,m,0,n,m);
  GiD = gsl_matrix_submatrix (G_inv,m,m,n,n);

  // block A and Schur complement of A
  gsl_matrix_memcpy (block_A, &GA.matrix);
  gsl_matrix_memcpy (Schur_A, &GD.matrix);

  // inversion of block A
  // = block A of inverse matrix (temporary)
  block_inv_4 (block_A, &GiA.matrix);

  // Ai_B (mxn) = A_inv (mxm) B (mxn)
  // C_Ai (nxm) = C (nxm) A_inv (mxm)
  gsl_matrix *Ai_B = gsl_matrix_alloc (m,n);
  gsl_matrix *C_Ai = gsl_matrix_alloc (n,m);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GiA.matrix,&GB.matrix,0,Ai_B);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GC.matrix,&GiA.matrix,0,C_Ai);

  // Shcur complement of block A (nxn) = C (nxm) A_inv (mxm) B (mxn)
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,C_Ai,&GB.matrix,1,Schur_A);

  // inversion of Shcur complement of block A (nxn) 
  // = block D of inverse matrix 
  block_inv_4 (Schur_A, &GiD.matrix);

  // block C of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,&GiD.matrix,C_Ai,0,&GiC.matrix);
  // block B of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiD.matrix,0,&GiB.matrix);
  // block A of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiC.matrix,1,&GiA.matrix);

  gsl_matrix_free (block_A);
  gsl_matrix_free (Schur_A);
  gsl_matrix_free (C_Ai);
  gsl_matrix_free (Ai_B);
}

void block_inv_4 (gsl_matrix *G, gsl_matrix *G_inv)
{
  size_t m = (size_t)(G->size1*0.5); 
  size_t n = G->size1 - m;
  gsl_permutation *qm = gsl_permutation_alloc (m);
  gsl_permutation *qn = gsl_permutation_alloc (n);
  gsl_matrix *block_A = gsl_matrix_alloc (m,m); // block A 
  gsl_matrix *Schur_A = gsl_matrix_alloc (n,n); // Schur complement of A
  gsl_matrix_view GA,GB,GC,GD,GiA,GiB,GiC,GiD;
  int signum;

  GA  = gsl_matrix_submatrix (G,0,0,m,m);
  GB  = gsl_matrix_submatrix (G,0,m,m,n);
  GC  = gsl_matrix_submatrix (G,m,0,n,m);
  GD  = gsl_matrix_submatrix (G,m,m,n,n);
  GiA = gsl_matrix_submatrix (G_inv,0,0,m,m);
  GiB = gsl_matrix_submatrix (G_inv,0,m,m,n);
  GiC = gsl_matrix_submatrix (G_inv,m,0,n,m);
  GiD = gsl_matrix_submatrix (G_inv,m,m,n,n);

  gsl_permutation_init (qm);
  gsl_permutation_init (qn);

  // block A and Schur complement of A
  gsl_matrix_memcpy (block_A, &GA.matrix);
  gsl_matrix_memcpy (Schur_A, &GD.matrix);

  // inversion of block A
  // = block A of inverse matrix (temporary)
  gsl_linalg_LU_decomp (block_A,qm,&signum);
  gsl_linalg_LU_invert (block_A,qm,&GiA.matrix);

  // Ai_B (mxn) = A_inv (mxm) B (mxn)
  // C_Ai (nxm) = C (nxm) A_inv (mxm)
  gsl_matrix *Ai_B = gsl_matrix_alloc (m,n);
  gsl_matrix *C_Ai = gsl_matrix_alloc (n,m);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GiA.matrix,&GB.matrix,0,Ai_B);
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&GC.matrix,&GiA.matrix,0,C_Ai);

  // Shcur complement of block A (nxn) = C (nxm) A_inv (mxm) B (mxn)
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,C_Ai,&GB.matrix,1,Schur_A);

  // inversion of Shcur complement of block A (nxn) 
  // = block D of inverse matrix 
  gsl_linalg_LU_decomp (Schur_A,qn,&signum);
  gsl_linalg_LU_invert (Schur_A,qn,&GiD.matrix);

  // block C of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,&GiD.matrix,C_Ai,0,&GiC.matrix);
  // block B of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiD.matrix,0,&GiB.matrix);
  // block A of inverse matrix
  gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,-1,Ai_B,&GiC.matrix,1,&GiA.matrix);

  gsl_matrix_free (block_A);
  gsl_matrix_free (Schur_A);
  gsl_matrix_free (C_Ai);
  gsl_matrix_free (Ai_B);
  gsl_permutation_free (qm);
  gsl_permutation_free (qn);
}
