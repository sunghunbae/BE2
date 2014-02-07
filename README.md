What is BE2?
============

BE2 is an ensemble approach to the boundary element method (BEM).
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
By default, an executable binary files (```be2```,```msms2gts```,```eg```) and 
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

MSMS accepts atomic coordinates as xyzr format 
which can be generated from a PDB file by an Awk script ```pdb_to_xyzr```.
The ```pdb_to_xyzr``` is a part of MSMS binary distribution and 
requires another file ```atmtypenumbers```. 
Note that ```atmtypenumbers``` should be in the current working directory.

```
$ chmod +x msms.i86Linux2.2.6.1
$ chmod +x pdb_to_xyzr
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
of 3*face columns and 3*face rows. BE2 calculation needs matrix inversion
and it may take a while for a big matrix. So, it would be sensible to reduce
the matrix size as long as the accuracy is not compromised severely.

The number of faces or elements depends on the size and shape complexitiy of 
the molecule of interest, but normally 800-2400 faces are good enough. 
However, best practice is to generate a series of surfaces with 
varying number of faces and run each of them, and extrapolate the results
to an infinite number of faces.

msms2gts serves this purpose. It reads MSMS output (.vert and .face files)
and coarsen it to a desired number of faces.

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

How to run for disordered ensemble?
===================================

### Generate ensemble

### Generate triangular surface

### Calculate diffusion tensors
