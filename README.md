What is BE2?
============

**BE2** is an ensemble approach to the boundary element method (BEM).
Using two layers of molecular surfaces whose correlated velocities decay exponentially with distance, 
it can accurately predicts molecular tumbling time for unstructed or disordered proteins.

### Background

For well-structured, rigid proteins, the prediction of rotational tumbling time using atomic coordinates is reasonably accurate, but is inaccurate for proteins with long unstructured sequences. Under physiological conditions, many proteins contain long disordered segments that play important regulatory roles in fundamental biological events including signal transduction and molecular recognition. 

Reliable prediction of tau(c) will help to detect intra- and intermolecular interactions and conformational switches between more ordered and less ordered states of the disordered segments. The method has been extensively validated using 12 reference proteins with 14 to 103 disordered residues at the N- and/or C-terminus and has been successfully employed to explain a set of published results on a system that incorporates a conformational switch.

### Citation

Please cite the following paper if you use BE2.

- Prediction of the Rotational Tumbling Time for Proteins with Disordered Segments.
  Sung-Hun Bae, H. Jane Dyson, and Peter E. Wright (2009)
  J. Am. Chem. Soc. 131 (19), 6814-6821.
  http://pubs.acs.org/doi/abs/10.1021/ja809687r

How to install?
===============

### Dependencies

You need the following packages to install BE2.

- CMake
- GNU Scientific Library (GSL) development
- GNU Triangulated Surface Library (GTS) development
- GLIB2.0 development

In Debian and related linux distributions (Ubuntu, Mint, etc):

```
$ sudo apt-get install cmake libgsl0-dev libgts-dev libglib2.0-dev
```

### Make
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make install
```
By default, an executable binary files (**```be2```**,**```msms2gts```**,**```eg```**) and 
data file directory (```eglib/```) will be installed in ~/bin directory. 
You may change the destination directory defined the CMakeLists.txt file :
```
SET(CMAKE_INSTALL_PREFIX ~/bin )
```

Please make sure that files and data are accessible from your working directory 
by adding ```~/bin```, ```~/bin/eglib/` or your destination directory in the 'path'.


How to run for rigid proteins?
==============================

### Generate triangular surface

Molecular surface represented by triangular surfaces are 
generated from atomic coordinates by MSMS program written by Michel F. Sanner. 
Binaries can be downloaded from http://mgltools.scripps.edu/downloads#msms

Please note that MSMS and pdb_to_xyzr should be given executable permission:

```
$ chmod +x msms.i86Linux2.2.6.1
$ chmod +x pdb_to_xyzr
```

MSMS accepts atomic coordinates as xyzr format 
which can be generated from a PDB file by an Awk script ```pdb_to_xyzr```.
The ```pdb_to_xyzr``` is a part of MSMS binary distribution and 
requires another file ```atmtypenumbers```. 
Note that ```atmtypenumbers``` should be in the current working directory.

```
$ pdb_to_xyzr 1UBQ.pdb | awk '{print $1,$2,$3,$4+1.1}' > 1UBQ.xyzr
$ msms.i86Linux2.2.6.1 -density 1 -if 1UBQ.xyzr -of 1UBQ >& msms.log
$ ls 1UBQ.*
  1UBQ.face 1UBQ.pdb 1UBQ.vert 1UBQ.xyzr
```

Here, we added 1.1 angstrom to the van der Waals radii to account for the 
hydration shell, which is commonly adopted practice in the shaped based calculation 
of hydrodynamic properties to make the calculation agree to the experimental
measurements.

BE2 accepts MSMS, GTS, and OFF format files. 
You can use the MSMS outputs (.vert and .face files) directly 
or convert them into equivalent or coarsened surfaces in the GTS format.
In the above ubiquitin example (PDBId: 1UBQ), MSMS would generate 
7970 triangular faces and 3987 vertices.

In BE2, these triangular faces or boundary elements are handled by a matrix 
of [3*face,3*face], thus it may take a long time to invert a large matrix.
For the above ubiquitin example, inversion of a 23970 by 23970 matrix derived
directly from MSMS outputs is not trivial and will take a very long time.
So, it would be sensible to minimize the matrix size as long as the 
accuracy is not severely compromised.

The appropriate number of faces or elements depends on the size and shape complexitiy of 
the molecule of interest, but normally 800-2400 faces are good enough for most protein
structures. However, best practice is to generate a series of surfaces with 
varying number of faces and run each of them, and extrapolate the results
to an infinite number of faces.

**msms2gts** was written to serve this purpose. 
It reads the MSMS outputs (.vert and .face files)
and coarsens them and saves to a GTS output with a desired number of faces.

```
$ msms2gts 1UBQ 800
$ ls 1UBQ.*
  1UBQ.0800.gts 1UBQ.face 1UBQ.pdb 1UBQ.vert 1UBQ.xyzr
```

### Calculate diffusion tensors

BE2 reads MSMS, GTS, or OFF format files and calculates translational and diffusion tensors.

```
$ be2 -gts 1UBQ > 1UBQ.0800.out
```

The output will be like:
```
Temperature  : 293.15 K
Viscosity    : 1.002 cP
GTS_surface  : 1UBQ.0800.gts (Vertex: 402 Edge: 1200 Face: 800)
Dimension    : X   36.4160  Y   36.9787 Z   42.1620
Calculating G matrix (2400 x 2400) ...... time 00:00:21
Inverting   G matrix (2400 x 2400) ...... time 00:01:19
Surface_Area (  0.534108 ...  12.356632) Sum= 4142.26 A^2

Center_of_Diffusion 30.7833 29.5617 16.6655
Dtt Eigenvalues (10^-7 cm^2 s^-1) | Eigenvectors
Dtt 1     12.281 |  8.377e-01 -2.428e-01 -4.891e-01
Dtt 2     12.527 | -4.117e-04 -8.960e-01  4.441e-01
Dtt 3     13.394 |  5.461e-01  3.719e-01  7.507e-01
Drr Eigenvalues (10^7 s^-1) | Eigenvectors
Drr 1      2.773 |  8.158e-01 -4.130e-01 -4.049e-01
Drr 2      2.890 | -1.374e-01 -8.184e-01  5.580e-01
Drr 3      4.009 |  5.618e-01  3.995e-01  7.244e-01
# rigid 1/(6Diso)   5.170 ns
```

### Accuracy

| PDB  | BE2 | Exp | PDB | BE2  | Exp  | PDB  |  BE2 |  Exp | PDB |  BE2  |  Exp  |
|------|-----|-----|-----|------|------|------|------|------|-----|-------|-------|
| 1ZNF | 2.0 | 2.4 |2CDS | 7.9  | 8.3  | 6I1B | 11.1 | 12.4 |1NWK | 25.2  | 24.5  |
| 5PTI | 4.3 | 4.0 |2CDS | 7.9  | 9.8  | 1LKI | 11.8 | 14.9 |1GKB | 31.7  | 32.7  |
| 1PIT | 4.4 | 4.4 |1HFX | 8.1  | 8.9  | 1TPO | 12.2 | 14.4 |1HHO | 33.5  | 29.8  |
| 2BCA | 4.5 | 5.1 |1PNE | 8.4  | 6.7  | 1SVN | 12.7 | 12.4 |2HHB | 34.7  | 32.7  |
| 1CLB | 4.9 | 4.9 |7RSA | 8.7  | 8.3  | 1BVG | 12.7 | 13.2 |1AO6 | 48.9  | 47.6  |
| 1UBQ | 5.2 | 5.4 |1E8L | 8.8  | 8.3  | 2CGA | 14.0 | 13.9 |1ALK | 53.1  | 53.8  |
| 1BTA | 5.6 | 7.4 |1MBO | 10.1 | 10.0 | 2CAB | 15.1 | 15.4 |2CTV | 60.8  | 59.5  |
| 1EGL | 5.7 | 6.2 |3BLG | 10.1 | 10.4 | 4PEP | 19.7 | 17.8 |6LDH | 86.0  | 83.3  |
| 1HRC | 6.7 | 6.9 |1AVU | 10.8 | 12.6 | 1BEB | 22.9 | 21.9 |1GPB | 125.1 | 128.2 |
| 2CDS | 7.9 | 7.6 |1WRT | 23.5 | 23.1 |

How to run for disordered proteins?
====================================

The disordered or unstructured proteins have one or more part(s) or domain(s)
that remain unstable and are constantly changing their conformations 
so that overall molecular shape fluctuates over time.

BE2 approach focuses on one rigid domain at a time.
However, it does not mean that BE2 is limited to proteins with only one rigid domain. 
For proteins with multiple rigid domains linked by some disordered linker(s), the 
diffusion tensor of each rigid domain can be separately calculated.
The rigid domain can be any size starting from just one residue. 

Since deposited PDB coordinates do not generally include the disordered 
flexible domain(s), flexible domain(s) should be modeled while retaining the 
structured domain(s). There are many ways to generate ensemble structures 
including the molecular modeling toolkit (MMTK).

A set of structures or ensemble describes the disordered proteins.
You may use your own ensemble structures or use **```eg```** to generate one.

**```eg```** generates ensemble structures by rotating backbone and side chain 
dihedral angles of the modeled template structure.
From the junction between the rigid and disordered domain(s) 
toward the N- or C-terminus, the molecular coordinates were
rotated according to a pair of backbone dihedral angles (phi-psi) and
side-chain dihedral angles (chi1-chi4) which were randomly selected
from an amino acid specific dihedral angle library. The phi-psi angle
library was built from 500 low-homology and high-resolution X-ray
structures (resolutions < 1.8 angstrom) excluding all residues in
alpha-helices, beta-sheets, and turns determined by the DSSP.
The chi angles library was adopted from a published rotamer library.
Residues immediately preceding proline were treated as an 
additional amino acid type, due to the restricted local conformation.
A generated ensemble structure is accepted only if
the number of van der Waals clash and van der Waals energy is less than the
given maxima. van der Waals energy is defined as:

> | Let rij and vdwij be the distance and sum of van der Waals radii of atoms i and j,  
| Evdw = 0 if vdwij < rij 
| Evdw = -57.273*(1.0-rij/(0.85*vdwio)) if  0.7*vdwij <= rij <= vdwij  
| Evdw = 10 if rij < 0.7*vdwij  

### Generate ensemble


### Generate triangular surface

### Calculate diffusion tensors
