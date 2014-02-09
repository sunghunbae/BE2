#include "rotate.h"

using namespace std;

int read_rot_list (const char *rot_file, vector < struct rot > &rlist)
{
    vector <string> c;
    string buffer;
    struct rot i;

    ifstream rfile (rot_file);

    if (rfile.is_open()) {
	while (!rfile.eof()) {
	    getline (rfile, buffer);
	    if (parse(buffer.c_str(), "\t ", c) != 0) {
		if (c.size() == 3) {
		    i.library_search = true;
		    sscanf (c[0].c_str(), "%c", &i.chain);
		    sscanf (c[1].c_str(), "%d", &i.resid);
		    sscanf (c[2].c_str(), "%c", &i.rotation);
		    rlist.push_back (i);
		    }
		if (c.size() == 10) {
		    i.library_search = false;
		    sscanf (c[0].c_str(), "%c", &i.chain);
		    sscanf (c[1].c_str(), "%d", &i.resid);
		    sscanf (c[2].c_str(), "%c", &i.rotation);
		    i.seq = c[3];
		    i.dihed[0] = c[4];
		    i.dihed[1] = c[5];
		    i.dihed[2] = c[6];
		    i.dihed[3] = c[7];
		    i.dihed[4] = c[8];
		    i.dihed[5] = c[9];
		    rlist.push_back (i);
		    } // read user-defined dihedral angles
		} // if not comment
	    } // until eof
	} // if rfile is opened
    else {
	cout << "cannot open " << rot_file << "\n";
	exit(-1);
	}
    return rlist.size();
}

void write_rot_list (const char *prefix, const vector < struct rot > & rlist,
    const int idx, const int num_clash, const double Evdw)
{
    char filename [256];
    sprintf (filename, "%s_%04d.rot",prefix,idx);
    FILE * fileO = fopen (filename, "w");
    fprintf(fileO,"# clash %d Evdw %.2f\n",num_clash,Evdw);
    for (int r=0; r<rlist.size();r++)
	fprintf (fileO,"%c %d %c %s %s %s %s %s %s %s\n",
	    rlist[r].chain, 
	    rlist[r].resid, 
	    rlist[r].rotation,
	    rlist[r].seq.c_str(),
	    rlist[r].dihed[0].c_str(), 
	    rlist[r].dihed[1].c_str(), 
	    rlist[r].dihed[2].c_str(),
	    rlist[r].dihed[3].c_str(), 
	    rlist[r].dihed[4].c_str(), 
	    rlist[r].dihed[5].c_str());
    fclose (fileO);
}

void rotate (vector <ATOM> &coor, vector <struct rot> &rlist, 
    vector <struct lib> &fyc)
{
    extern gsl_rng * rng;
    int r,i,j,v,w;

    for (r=0; r< rlist.size(); r++) {

	if (rlist[r].library_search) {
    
	    string next_seq;

	    get_sequence (coor, rlist[r].chain, rlist[r].resid, rlist[r].seq);
	    get_sequence (coor, rlist[r].chain, rlist[r].resid+1, next_seq);

	    // select phi/psi dihedral angles randomly from library
	    for (i=0;i<fyc.size();i++) {
		if (rlist[r].seq == fyc[i].seq) {
		    if (next_seq == "PRO") {
			v = (int) gsl_rng_uniform_int (rng,fyc[i].fP.size());
			rlist[r].dihed[0] = fyc[i].fP[v];
			rlist[r].dihed[1] = fyc[i].yP[v];
			}
		    else {
			v = (int) gsl_rng_uniform_int (rng,fyc[i].f.size());
			rlist[r].dihed[0] = fyc[i].f[v];
			rlist[r].dihed[1] = fyc[i].y[v];
			}
		    break;
		    } // if
		} // for

	    // select chi dihedral angles randomly from library
	    if (fyc[i].c1.size() > 0 ) {
		w = (int) gsl_rng_uniform_int (rng,fyc[i].c1.size());
		rlist[r].dihed[2] = fyc[i].c1[w];
		rlist[r].dihed[3] = fyc[i].c2[w];
		rlist[r].dihed[4] = fyc[i].c3[w];
		rlist[r].dihed[5] = fyc[i].c4[w];
		}
	    else {
		rlist[r].dihed[2] = "-";
		rlist[r].dihed[3] = "-";
		rlist[r].dihed[4] = "-";
		rlist[r].dihed[5] = "-";
		}
	    
	    } // library search

	// rotate coordinates according to selected dihedral angles
	dihedrot (coor, rlist[r]);

	} // r
}

void dihedrot (vector <ATOM> & coor, const struct rot &R)
{
    Topology TOPO;
    Vector _n,_nn,_cp,_c,_a,_b,_g,_d,_e,_z;
    vector <string> p;
    double diff;
    int i;

    // rotate phi
    if (R.dihed[0] != "-") {
	get_vector (coor, R.chain, R.resid-1, "C", _cp);
	get_vector (coor, R.chain, R.resid, "N", _n);
	get_vector (coor, R.chain, R.resid, "CA", _a);
	get_vector (coor, R.chain, R.resid, "C", _c);
	diff = Dihedral_ (&_cp,&_n,&_a,&_c) - strtod (R.dihed[0].c_str(),NULL);
	for (i=0;i<coor.size();i++) {
	    if (coor[i].chainId == R.chain) {
		/* rotate N-terminal part */
		if (R.rotation == 'N' && (
		    (coor[i].resSeq < R.resid) || 
		    (coor[i].resSeq == R.resid && coor[i].name_ == "HN")))
		    RotateAtom_ (coor[i],&_n,&_a,diff);
		/* rotate C-terminal part */
		if (R.rotation == 'C' && (
		    (coor[i].resSeq > R.resid) || 
		    (coor[i].resSeq == R.resid && coor[i].name_ != "HN")))
		    RotateAtom_ (coor[i],&_a,&_n,diff);
		} // within the same chain
	    } 
	}

    // rotate psi
    if (R.dihed[1] != "-") {
	get_vector(coor, R.chain, R.resid, "N", _n);
	get_vector(coor, R.chain, R.resid, "CA",_a);
	get_vector(coor, R.chain, R.resid, "C", _c);
	get_vector(coor, R.chain, R.resid+1, "N", _nn);
	diff = Dihedral_ (&_n,&_a,&_c,&_nn) - strtod(R.dihed[1].c_str(),NULL);
	for (i=0;i<coor.size();i++) {
	    if (coor[i].chainId == R.chain) {
		/* rotate N-terminal part */
		if (R.rotation == 'N' && (
		    (coor[i].resSeq < R.resid) || 
		    (coor[i].resSeq == R.resid && coor[i].name_ != "O")))
		    RotateAtom_ (coor[i],&_a,&_c,diff);
		/* rotate C-terminal part */
		if (R.rotation == 'C' && (
		    (coor[i].resSeq > R.resid) || 
		    (coor[i].resSeq == R.resid && coor[i].name_ == "O")))
		    RotateAtom_ (coor[i],&_c,&_a,diff);
		} // within the same chain
	    }
	}

    // Chi 1
    if(R.dihed[2] != "-" && TOPO.chi(1, R.seq)) {
	TOPO.chi(1, R.seq, p);
	get_vector (coor, R.chain, R.resid, p[0],_n);
	get_vector (coor, R.chain, R.resid, p[1],_a);
	get_vector (coor, R.chain, R.resid, p[2],_b);
	get_vector (coor, R.chain, R.resid, p[3],_g);
	diff = Dihedral_ (&_n,&_a,&_b,&_g) - strtod(R.dihed[2].c_str(),NULL);
	for (i=0; i < coor.size(); i++) {
	    if (coor[i].chainId == R.chain) {
		if(coor[i].resSeq == R.resid && 
		    coor[i].name_ != "HN" && coor[i].name_ != "N" &&
		    coor[i].name_ != "C" && coor[i].name_ != "O" &&
		    coor[i].name_ != "CA" && coor[i].name_ != "HA")
		    RotateAtom_ (coor[i],&_b,&_a,diff);
		} // within the same chain
	    }
	}

    // Chi 2
    if(R.dihed[3] != "-" && TOPO.chi(2, R.seq)) {
	TOPO.chi(2,R.seq,p); 
	get_vector (coor, R.chain, R.resid, p[0],_a);
	get_vector (coor, R.chain, R.resid, p[1],_b);
	get_vector (coor, R.chain, R.resid, p[2],_g);
	get_vector (coor, R.chain, R.resid, p[3],_d);
	diff = Dihedral_ (&_a,&_b,&_g,&_d) - strtod(R.dihed[3].c_str(),NULL);
	for (i=0; i < coor.size(); i++) {
	    if (coor[i].chainId == R.chain) {
		if( coor[i].resSeq == R.resid && 
		    coor[i].name_ != "HN" && coor[i].name_ != "N" &&
		    coor[i].name_ != "C" && coor[i].name_ != "O" &&
		    coor[i].name_ != "CA" && coor[i].name_ != "HA" &&
		    coor[i].name_ != "CB" && coor[i].name_ != "HB1" &&
		    coor[i].name_ != "HB2" )
		    RotateAtom_ (coor[i],&_g,&_b,diff);
		}// within the same chain
	    }
	}

    // Chi 3
    if(R.dihed[4] != "-" && TOPO.chi(3,R.seq)) {
	TOPO.chi (3,R.seq,p); 
	get_vector (coor, R.chain, R.resid, p[0],_b);
	get_vector (coor, R.chain, R.resid, p[1],_g);
	get_vector (coor, R.chain, R.resid, p[2],_d);
	get_vector (coor, R.chain, R.resid, p[3],_e);
	diff = Dihedral_ (&_b,&_g,&_d,&_e) - strtod(R.dihed[4].c_str(),NULL);
	for (i=0; i < coor.size(); i++) {
	    if (coor[i].chainId == R.chain) {
		if(coor[i].resSeq == R.resid && 
		    coor[i].name_ != "HN" && coor[i].name_ != "N" &&
		    coor[i].name_ != "C" && coor[i].name_ != "O" &&
		    coor[i].name_ != "CA" && coor[i].name_ != "HA" &&
		    coor[i].name_ != "CB" && coor[i].name_ != "HB1" &&
		    coor[i].name_ != "HB2" && coor[i].name_ != "CG" &&
		    coor[i].name_ != "HG1" && coor[i].name_ != "HG2")
		    RotateAtom_ (coor[i],&_d,&_g,diff);
		}// within the same chain
	    }
	}

    // Chi 4
    if(R.dihed[5] != "-" && TOPO.chi(4, R.seq)) {
	TOPO.chi (4,R.seq,p);
	get_vector (coor, R.chain, R.resid, p[0], _g);
	get_vector (coor, R.chain, R.resid, p[1], _d);
	get_vector (coor, R.chain, R.resid, p[2], _e);
	get_vector (coor, R.chain, R.resid, p[3], _z);
	diff = Dihedral_ (&_g,&_d,&_e,&_z) - strtod(R.dihed[5].c_str(),NULL);
	for (i=0; i< coor.size(); i++) {
	    if (coor[i].chainId == R.chain) {
		if(coor[i].resSeq == R.resid && 
		    coor[i].name_ != "HN" && coor[i].name_ != "N" &&
		    coor[i].name_ != "C" && coor[i].name_ != "O" &&
		    coor[i].name_ != "CA" && coor[i].name_ != "HA" &&
		    coor[i].name_ != "CB" && coor[i].name_ != "HB1" &&
		    coor[i].name_ != "HB2" && coor[i].name_ != "CG" &&
		    coor[i].name_ != "HG1" && coor[i].name_ != "HG2" && 
		    coor[i].name_ != "CD" && coor[i].name_ != "HD1" &&
		    coor[i].name_ != "HD2" )
		    RotateAtom_ (coor[i],&_e,&_d,diff);
		}// within the same chain
	    }
	}
}


void get_sequence (const vector <ATOM> &coor,
  const char chainId, const size_t resSeq, string &resName)
{
  resName = "";
  for (size_t i=0; i< coor.size(); i++) {
    if (coor[i].chainId == chainId && coor[i].resSeq == resSeq) {
      resName = coor[i].resName;
      break;
      }
    }
}

void get_vector (const vector <ATOM> &coor,
  char chainId, size_t resSeq, string p, Vector &a)
{
    for (size_t i=0; i< coor.size(); i++) {
	if (coor[i].chainId == chainId && 
	    coor[i].resSeq == resSeq && coor[i].name_ == p) {
	    a.x = coor[i].x;
	    a.y = coor[i].y;
	    a.z = coor[i].z;
	    break;
	    }
	}
}

int RotateAtom_ (ATOM & atom, Vector *v1, Vector *v2, double angle)
{
  Vector rot_axis_vectorS;
  Vector radius_vectorS;
  Vector parallel_vectorS;
  Vector perpendicular_vectorS;
  Vector unit_vector1S, unit_vector2S;
  double  abs_value;
  double  perpendicular_part;
  double  reciprocal_abs_value;
  double  cos_angle, sin_angle;
  double  p1_new, p2_new;

  /* Prepare the vector which defines the rotation axis: */
  rot_axis_vectorS.x = v2->x - v1->x;
  rot_axis_vectorS.y = v2->y - v1->y;
  rot_axis_vectorS.z = v2->z - v1->z;

  /* Check the absolute value of the rotation axis vector: */
  abs_value = AbsoluteValue_ (&rot_axis_vectorS);
  if (abs_value <= 0.0) return -1;

  /* Prepare the radius vector of the given atom: */
  radius_vectorS.x = atom.x - v1->x;
  radius_vectorS.y = atom.y - v1->y;
  radius_vectorS.z = atom.z - v1->z;

  /* Find the part of the radius vector which is parallel to the rotation */
  /* axis.  Note that  we need the vector,  not just  the absolute value! */
  ParallelPart_ (&parallel_vectorS, &rot_axis_vectorS, &radius_vectorS);

  /* Find the part of the radius vector which */
  /* is perpendicular  to the  rotation axis: */
  perpendicular_vectorS.x = radius_vectorS.x - parallel_vectorS.x;
  perpendicular_vectorS.y = radius_vectorS.y - parallel_vectorS.y;
  perpendicular_vectorS.z = radius_vectorS.z - parallel_vectorS.z;

  /* Prepare and check the absolute value of the perpendicular part: */
  perpendicular_part = AbsoluteValue_ (&perpendicular_vectorS);
  if (perpendicular_part <= 0.0) return -2;

  /* Prepare the first unit vector, required for rotation: */
  reciprocal_abs_value = 1.0 / perpendicular_part;
  unit_vector1S.x = reciprocal_abs_value * perpendicular_vectorS.x;
  unit_vector1S.y = reciprocal_abs_value * perpendicular_vectorS.y;
  unit_vector1S.z = reciprocal_abs_value * perpendicular_vectorS.z;

  /* Prepare and check the second unit vector: */
  VectorProduct_ (&unit_vector2S, &rot_axis_vectorS, &unit_vector1S);
  abs_value = AbsoluteValue_ (&unit_vector2S);
  if (abs_value <= 0.0) return -3;
  reciprocal_abs_value = 1.0 / abs_value;
  unit_vector2S.x = reciprocal_abs_value * unit_vector2S.x;
  unit_vector2S.y = reciprocal_abs_value * unit_vector2S.y;
  unit_vector2S.z = reciprocal_abs_value * unit_vector2S.z;

  /* Rotate the perpendicular vector: */
  cos_angle = cos (angle*M_PI/180.0);
  sin_angle = sin (angle*M_PI/180.0);
  p1_new = perpendicular_part * cos_angle;
  p2_new = perpendicular_part * sin_angle;
  perpendicular_vectorS.x = p1_new * unit_vector1S.x + p2_new * unit_vector2S.x;
  perpendicular_vectorS.y = p1_new * unit_vector1S.y + p2_new * unit_vector2S.y;
  perpendicular_vectorS.z = p1_new * unit_vector1S.z + p2_new * unit_vector2S.z;

  /* Update the radius vector: */
  radius_vectorS.x = parallel_vectorS.x + perpendicular_vectorS.x;
  radius_vectorS.y = parallel_vectorS.y + perpendicular_vectorS.y;
  radius_vectorS.z = parallel_vectorS.z + perpendicular_vectorS.z;

  /* Update the atomic coordinates: */
  atom.x = radius_vectorS.x + v1->x;
  atom.y = radius_vectorS.y + v1->y;
  atom.z = radius_vectorS.z + v1->z;
  
  /* Return positive value (success indicator): */
  return 1;
}


double Dihedral_ (Vector *a, Vector *b, Vector *c, Vector *d)
{
  Vector ba, bc, cb, cd, u1S, u2S, v1S, v2S;
  double  denom, ratio, alpha, dihe;

  ba.x = a->x - b->x;
  ba.y = a->y - b->y;
  ba.z = a->z - b->z;
  bc.x = c->x - b->x;
  bc.y = c->y - b->y;
  bc.z = c->z - b->z;

  cb.x = b->x - c->x;
  cb.y = b->y - c->y;
  cb.z = b->z - c->z;
  cd.x = d->x - c->x;
  cd.y = d->y - c->y;
  cd.z = d->z - c->z;

  /* Two vectors perpendicular to bc vector,  mutually orthogonal, the */
  /* second in the plane defined by N_previousC_vectorS and N_CA_vectorS: */
  VectorProduct_ (&u1S, &ba, &bc);
  VectorProduct_ (&u2S, &u1S, &bc);

  /* Two vectors perpendicular to  CA_N_vectorS,  mutually orthogonal, */
  /* the second in the plane defined by CA_N_vectorS and CA_C_vectorS: */
  VectorProduct_ (&v1S, &cb, &cd);
  VectorProduct_ (&v2S, &cb, &v1S);

  /* Calculate the angle alpha, which will be used to calculate phi: */
  /* Avoid division by zero: */
  denom = AbsoluteValue_ (&u1S) * AbsoluteValue_ (&v1S);
  if (denom == 0.0) return -999;
  /* Use the scalar product to calculate the cosine of the angle: */
  ratio = ScalarProduct_ (&u1S, &v1S) / denom;
  /* Arc cosine is very sensitive to floating point errors: */
  if (ratio <= -1.0) alpha = 3.1415927;
  else if (ratio >= 1.0) alpha = 0.0;
  else alpha = acos (ratio);
  /* There are two possible solutions; the right one is resolved here: */
  if (ScalarProduct_ (&v2S, &u1S) >= 0) dihe = alpha;
  else dihe = -alpha;
  /* Return the angle (in degree): */
  return dihe*180.0/M_PI;
}


void VectorProduct_ (Vector *v3,Vector *v1,Vector *v2)
{
  v3->x = v1->y * v2->z - v1->z * v2->y;
  v3->y = v1->z * v2->x - v1->x * v2->z;
  v3->z = v1->x * v2->y - v1->y * v2->x;
}

double ScalarProduct_ (Vector *v1, Vector *v2) 
{
  return (v1->x * v2->x + v1->y * v2->y + v1->z * v2->z); 
}

double AbsoluteValue_ (Vector *v) 
{
  return sqrt (v->x * v->x + v->y * v->y + v->z * v->z); 
}

int ParallelPart_ ( Vector *parallel_v, Vector *reference_v, Vector *input_v) 
{
  double  abs_value_squared;
  double  reciprocal_denominator;
  double  scalar_product;
  double  scale_factor;
  abs_value_squared = ScalarProduct_ (reference_v, reference_v);
  if (abs_value_squared == 0.0) return -1;
  reciprocal_denominator = 1.0 / abs_value_squared;
  scalar_product = ScalarProduct_ (input_v, reference_v);
  scale_factor = scalar_product * reciprocal_denominator;
  parallel_v->x = scale_factor * reference_v->x;
  parallel_v->y = scale_factor * reference_v->y;
  parallel_v->z = scale_factor * reference_v->z;
  return 1; 
}
