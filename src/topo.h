#ifndef __TOPO__
#define __TOPO__

#include <string>
#include <vector>
#include <iostream>

using namespace std;

typedef struct {
    string seq;
    vector < vector < string > > dihedral;
    } aminoacid;

class Topology {
    private:
	vector < aminoacid > topo;
	void add_chi (const string seq, 
		const string, const string, const string, const string);
    public:
	Topology();
	~Topology ();
	bool chi (const int, const string seq);
	void chi (const int, const string seq, vector < string > &);
	void view_chi ();
};

#endif
