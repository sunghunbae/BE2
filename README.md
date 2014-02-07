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
-- Ubuntu: sudo apt-get install libgsl0-dev
- GNU Triangulated Surface Library (GTS) development
- GLIB2.0 development


### Make
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make install
```

How to run?
===========

### Generate triangular surfaces from PDB

We need MSMS program written by Michel F. Sanner to compute molecular surfaces 
from atomic coordinates. Binaries can be downloaded from http://mgltools.scripps.edu/downloads#msms

```
$ pdb_to_xyzr a.pdb | awk '{print $1,$2,$3,$4+1.1}' > a.xyzr
```
Here, we assumed 1.1 angstrom hydration shell 
which is commonly adopted in the shaped based calculation of hydrodynamic properties.
MSMS then generates triangular surfaces from the xyzr format.

```
$ msms.i86Linux2.2.6.1 -density 1 -if a.xyzr -of a >& msms.log
```

BE2 accepts MSMS, GTS, and OFF format files. You may directly use the MSMS outputs 
(.vert and .face files) or convert them into equivalent or coarsened surfaces in the GTS format.

```
$ msms2gts a 1600
```

BE2 performs matrix inversion and it may take a long time if the number of faces
are more than necessarily large. It depends on the size and shape of the molecule of interest,
but normally 800-2400 faces are good enough.

### Calculate diffusion tensors

```
$ be2 -gts a > a.out
```


