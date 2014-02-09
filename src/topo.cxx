#include "topo.h"

using namespace std;

bool Topology::chi (const int idx, const string seq)
{
    int i;
    for (i=0; i < topo.size(); i++) {
	if (topo[i].seq == seq && topo[i].dihedral.size() >= idx &&
	    topo[i].dihedral[idx-1].size() == 4)
		return true;
	} 
    return false;
}

void Topology::chi (const int idx, const string seq, vector < string > &p)
{
    int i;
    for (i=0; i < topo.size(); i++) {
	if (topo[i].seq == seq && topo[i].dihedral.size() >= idx &&
	    topo[i].dihedral[idx-1].size() == 4) {
		p = topo[i].dihedral[idx-1];
		return;
		}
	} 
}

void Topology::add_chi (const string seq, const string p, 
    const string q, const string r, const string s)
{
    int i;
    bool found = false;
    vector < string > atoms;
    aminoacid a;

    atoms.empty();
    atoms.push_back (p);
    atoms.push_back (q);
    atoms.push_back (r);
    atoms.push_back (s);

    for (i = 0; i < topo.size(); i++)
	if (topo[i].seq == seq) {
	    found = true;
	    break;
	    }

    if (found) {
	topo[i].dihedral.push_back (atoms);
	}
    else {
	a.seq = seq;
	a.dihedral.push_back (atoms);
	topo.push_back (a);
	}
}

void Topology::view_chi ()
{
    int i,j,k;
    for (i=0;i<topo.size();i++) {
	for (j=0;j<topo[i].dihedral.size();j++) {
	    cout << "TOPOLOGY:"<<topo[i].seq<<"-chi"<<(j+1)<<" ";
	    for (k=0;k<topo[i].dihedral[j].size();k++) {
		cout << topo[i].dihedral[j][k] << " ";
		}
	    cout << endl;
	    }
	}
}

Topology::~Topology() { topo.clear(); }
Topology::Topology ()
{
    // chi1
    add_chi ("ARG", "N","CA","CB","CG");
    add_chi ("ASN", "N","CA","CB","CG");
    add_chi ("ASP", "N","CA","CB","CG");
    add_chi ("CYS", "N","CA","CB","SG");
    add_chi ("GLN", "N","CA","CB","CG");
    add_chi ("GLU", "N","CA","CB","CG");
    add_chi ("HIS", "N","CA","CB","CG");
    add_chi ("ILE", "N","CA","CB","CG1");
    add_chi ("LEU", "N","CA","CB","CG");
    add_chi ("LYS", "N","CA","CB","CG");
    add_chi ("MET", "N","CA","CB","CG");
    add_chi ("PHE", "N","CA","CB","CG");
    add_chi ("PRO", "N","CA","CB","CG");
    add_chi ("SER", "N","CA","CB","OG");
    add_chi ("THR", "N","CA","CB","OG1");
    add_chi ("TRP", "N","CA","CB","CG");
    add_chi ("TYR", "N","CA","CB","CG");
    add_chi ("VAL", "N","CA","CB","CG1");

    // chi2
    add_chi ("ARG","CA","CB","CG","CD");
    add_chi ("ASN","CA","CB","CG","OD1");
    add_chi ("ASP","CA","CB","CG","OD1");
    add_chi ("GLN","CA","CB","CG","CD");
    add_chi ("GLU","CA","CB","CG","CD");
    add_chi ("HIS","CA","CB","CG","ND1");
    add_chi ("ILE","CA","CB","CG1","CD1");
    add_chi ("LEU","CA","CB","CG","CD1");
    add_chi ("LYS","CA","CB","CG","CD");
    add_chi ("MET","CA","CB","CG","SD");
    add_chi ("PHE","CA","CB","CG","CD1");
    add_chi ("PRO","CA","CB","CG","CD");
    add_chi ("TRP","CA","CB","CG","CD1");
    add_chi ("TYR","CA","CB","CG","CD1");

    // chi3
    add_chi ("ARG","CB","CG","CD","NE");
    add_chi ("GLN","CB","CG","CD","OE1");
    add_chi ("GLU","CB","CG","CD","OE1");
    add_chi ("LYS","CB","CG","CD","CE");
    add_chi ("MET","CB","CG","SD","CE");

    // chi4
    add_chi ("ARG","CG","CD","NE","CZ");
    add_chi ("LYS","CG","CD","CE","NZ");

    // chi5
    add_chi ("ARG","CD","NE","CZ","NH1");
}

/* 
	DEFINITION OF DIHEDRAL ANGLES

+-----+-------------+--------------+--------------+-------------+--------------+
| res.|    chi1     |     chi2     |     chi3     |     chi4    |     chi5     |
|-----+-------------+--------------+--------------+-------------+--------------|
| ALA |             |              |              |             |              |
| ARG | N-CA-CB-CG  | CA-CB-CG-CD  | CB-CG-CD-NE  | CG-CD-NE-CZ | CD-NE-CZ-NH1 |
| ASN | N-CA-CB-CG  | CA-CB-CG-OD1 |              |             |              |
| ASP | N-CA-CB-CG  | CA-CB-CG-OD1 |              |             |              |
| CYS | N-CA-CB-SG  |              |              |             |              |
| GLN | N-CA-CB-CG  | CA-CB-CG-CD  | CB-CG-CD-OE1 |             |              |
| GLU | N-CA-CB-CG  | CA-CB-CG-CD  | CB-CG-CD-OE1 |             |              |
| GLY |             |              |              |             |              |
| HIS | N-CA-CB-CG  | CA-CB-CG-ND1 |              |             |              |
| ILE | N-CA-CB-CG1 | CA-CB-CG1-CD1|              |             |              |
| LEU | N-CA-CB-CG  | CA-CB-CG-CD1 |              |             |              |
| LYS | N-CA-CB-CG  | CA-CB-CG-CD  | CB-CG-CD-CE  | CG-CD-CE-NZ |              |
| MET | N-CA-CB-CG  | CA-CB-CG-SD  | CB-CG-SD-CE  |             |              |
| PHE | N-CA-CB-CG  | CA-CB-CG-CD1 |              |             |              |
| PRO | N-CA-CB-CG  | CA-CB-CG-CD  |              |             |              |
| SER | N-CA-CB-OG  |              |              |             |              |
| THR | N-CA-CB-OG1 |              |              |             |              |
| TRP | N-CA-CB-CG  | CA-CB-CG-CD1 |              |             |              |
| TYR | N-CA-CB-CG  | CA-CB-CG-CD1 |              |             |              |
| VAL | N-CA-CB-CG1 |              |              |             |              |
+-----+-------------+--------------+--------------+-------------+--------------+

*/
