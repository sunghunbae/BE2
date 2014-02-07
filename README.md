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

MSMS accepts xyzr format which can be generated from a PDB file.
```pdb_to_xyzr``` is a part of MSMS binary distribution and requires another
file ```atmtypenumbers```. Note that ```pdb_to_xyzr``` and 
```atmtypenumbers``` should be in the same directory and 
accessible from your working directory.

```
$ pdb_to_xyzr a.pdb | awk '{print $1,$2,$3,$4+1.1}' > a.xyzr
```

Here, we added 1.1 angstrom to the van der Waals radii to account for the 
hydration shell, which is commonly adopted practice in the shaped based calculation 
of hydrodynamic properties to make the calculation agree to the experimental
measurements.

Execution of the following command would result in two files, ```a.vert``` and ```a.face```.

```
$ msms.i86Linux2.2.6.1 -density 1 -if a.xyzr -of a >& msms.log
```

BE2 accepts MSMS, GTS, and OFF format files. 
You can use the MSMS outputs (.vert and .face files) directly 
or convert them into equivalent or coarsened surfaces in the GTS format.

BE2 performs matrix inversions and it may take a long time if the number of faces,
and therefore the size of matrix is more than necessarily large.
It depends on the size and shape of the molecule of interest,
but normally 800-2400 faces are good enough. However, best practice would be
generating a series of GTS format files with varying number of faces and
run each of them, and extrapolate to an infinite number of faces.

```
$ msms2gts a 1600
```

This will read a.vert and a.face and coarsen the triangular surfaces down to 1600 faces
and save to a.gts.

### Calculate diffusion tensors

```
$ be2 -gts a > a.out
```

BE2 reads MSMS, GTS, or OFF format files and calculates translational and diffusion tensors.
The output will be like:
```
GTS_surface  : 1NWK.1200.gts (Vertex: 600 Edge: 1800 Face: 1200)
Dimension    : X   75.9305  Y   53.2333 Z   59.3743
Temperature  : 293.15 K
Viscosity    : 1.002 cP
Calculating G matrix (3600 x 3600) ...... time 00:00:38
Inverting   G matrix (3600 x 3600) ...... time 00:10:01
Surface_Area (  0.774144 ...  35.811728) Sum= 13233.6 A^2
Center_Of_Diffusion (   8.6754,  -0.5887,  21.3645)

Eigenvalue (10^-7 cm^2 s^-1) & Eigenvector of Dtt
1      7.102 | -1.982e-01 -7.904e-01  5.796e-01
2      7.723 | -2.844e-01  6.123e-01  7.377e-01
3      7.936 | -9.380e-01 -1.860e-02 -3.461e-01

Eigenvalue (10^7 s^-1) & Eigenvector of Drr
1      0.564 | -2.109e-01 -7.994e-01  5.626e-01
2      0.664 | -3.129e-01  6.005e-01  7.359e-01
3      0.777 |  9.261e-01  2.088e-02  3.768e-01
# 1/(6Diso)  24.932 ns
```

How to run for disordered ensemble?
===================================

### Generate ensemble

### Generate triangular surface

### Calculate diffusion tensors
