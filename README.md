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

- Prediction of the Rotational Tumbling Time for Proteins with Disordered Segments
Sung-Hun Bae, H. Jane Dyson, and Peter E. Wright
Journal of the American Chemical Society 2009 131 (19), 6814-6821
http://pubs.acs.org/doi/abs/10.1021/ja809687r

How to install?
===============

### Dependencies

You need the following packages to install BE2.

- CMake
- GNU Scientific Library (GSL)
- GNU Triangulated Surface Library (GTS)

### Make
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make install
```
