What is BE2?
============

**BE2** is an ensemble approach to the boundary element method (BEM).
Using two layers of molecular surfaces whose correlated velocities decay exponentially with distance, 
it can accurately predict molecular tumbling time, tauc, for rigid and flexible or disordered proteins.

### Background

For well-structured rigid proteins, the prediction of rotational tumbling time using atomic coordinates is reasonably accurate, but is inaccurate for proteins with long unstructured sequences. Under physiological conditions, many proteins contain long disordered segments that play important regulatory roles in fundamental biological events including signal transduction and molecular recognition. 

Reliable prediction of tauc will help to detect intra- and intermolecular interactions and conformational switches between more ordered and less ordered states of the disordered segments. The method has been extensively validated using 12 reference proteins with 14 to 103 disordered residues at the N- and/or C-terminus and has been successfully employed to explain a set of published results on a system that incorporates a conformational switch.

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
a data file ```fycqr.lib``` will be installed at ```~/bin``` 
and ```~/bin/eglib``` directories, respectively. 
You may change the destination directory defined the CMakeLists.txt file :
```
SET(CMAKE_INSTALL_PREFIX ~/bin )
```

Please make sure that files are accessible from your working directory 
by adding ```~/bin```, ```~/bin/eglib/` or your destination directory in the ```path```.


How to run for the rigid proteins?
==================================

### Generate triangular surface

The triangular boundary elements representing molecular surface are 
generated from atomic coordinates by MSMS program written by Michel F. Sanner. 
MSMS binaries can be downloaded from http://mgltools.scripps.edu/downloads#msms

Please note that MSMS and pdb_to_xyzr should be given executable permission after download.

```
$ chmod +x msms.i86Linux2.2.6.1
$ chmod +x pdb_to_xyzr
```

MSMS accepts atomic coordinates as xyzr format 
which can be generated from a PDB file by an Awk script ```pdb_to_xyzr```.
The ```pdb_to_xyzr``` is a part of MSMS binary distribution and 
requires another file ```atmtypenumbers```. 
Note that ```pdb_to_xyzr``` looks for ```atmtypenumbers``` in the current directory(```.```).

```
$ pdb_to_xyzr 1UBQ.pdb | awk '{print $1,$2,$3,$4+1.1}' > 1UBQ.xyzr
$ msms.i86Linux2.2.6.1 -density 1 -if 1UBQ.xyzr -of 1UBQ >& msms.log
$ ls 1UBQ.*
  1UBQ.face 1UBQ.pdb 1UBQ.vert 1UBQ.xyzr
```

Here, we added 1.1 angstrom to the van der Waals radius of an atom to account for the 
arguable hydration shell, which is commonly assumed in the shaped based calculation 
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
So, it is sensible to minimize the matrix size as long as the 
accuracy is not severely compromised.

The appropriate number of faces or boundary elements depends on the size and shape 
complexitiy of the molecule of interest, 
but normally 800-2400 faces are good enough for most protein structures. 
However, best practice is to generate a series of surfaces with 
varying number of faces and run each of them, and extrapolate the results
to an infinite number of faces.

**msms2gts** serves this purpose. 
It reads and coarsens the MSMS outputs (.vert and .face files)
and saves to a GTS output with a desired number of faces.

```
$ msms2gts 1UBQ 800
$ ls 1UBQ.*
  1UBQ.0800.gts 1UBQ.face 1UBQ.pdb 1UBQ.vert 1UBQ.xyzr
```

### Calculate diffusion tensors

Running BE2 is very simple.

```
$ be2 -gts 1UBQ.0800
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

The accuracy of BE2 was verified with the rigid proteins and pure geometric objects 
(not shown here).

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

How to run for the disordered proteins?
=======================================

The disordered or flexible or unstructured proteins have 
one or more part(s) or domain(s)
that remain unstable and are constantly changing their conformations 
such that overall molecular shape can fluctuate over time.

Extending the classical boundary element method(BEM) applied to the rigid objects, 
BE2 approach focuses on one rigid domain while regarding all the other parts of 
a molecule as *environment* that is dragging motion in the solution.
However, it does not mean that BE2 is limited to proteins with only one rigid domain. 
For proteins with multiple rigid domains linked by some disordered or flexible linker(s), 
the diffusion tensors can be calculated individually for each domain.
It is also important to mention that the rigid domain can be any size. 
For example, a completely disordered chain of amino acid can be also investigated with BE2
by regarding one amino acid as the rigid domain at a time with the rest of amino acids being
the *environment*.

Generally, deposited PDB coordinates do not include the disordered or 
flexible domain(s). So, flexible domain(s) should be modeled while retaining the 
structured domain(s) to generate ensemble structures. 

You may use your own ensemble structures for BE2. 
Otherwise, you can generate a template structure of the whole sequence including
the disordered or flexible part(s) by using a modeling software such MMTK,
then run **```eg```**  to sample the conformational space.

### Generate ensemble

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
the number of van der Waals clash and van der Waals energy 
are below the given limits.

Please note that the rigid domain in the generated pdb files will have  
the same XYZ coordinates as the template pdb, which helps BE2 to relate
the rigid domain surface and the ensemble structure surfaces. 

```
$ eg
  Ensemble Generator (coil library/dihedral angle rotations)
  Usage: eg [options]
  options:
    -i #        read <i>nput PDB file
    -d # # #    rigid <d>omain [chainId begin end]
    -r #        read <r>otation list (default= rlist)
    -l #        import <l>ibrary (default= fycqr.lib)
    -o # # #    <o>utput [prefix begin end]
    -maxc #     <max>imum allowed <c>lashes (default= 0)
    -maxE #     <max>imum allowed VDW <e>nergy (default= 10000)
    -seed #     <seed> for random number generation
    -pdb        save output as <PDB> (.pdb)
    -rlist      save output as <rlist> (.rot)
    -rebuild    re<b>uild from rlist
    -quiet      be <quiet> do not print anything
```

The mouse prion protein(89-230) ensemble structures can be generated 
using the following *rlist* (rotation list) file.

```
$ cat rlist-prp
A 90 N
A 91 N
A 92 N
--- omitted ---
A 124 N
A 125 N
A 126 N
#
# folded-domain (127-224, coordinates fixed)
#
A 225 C
A 226 C
A 227 C
A 228 C
A 229 C
A 230 C

$ eg -i MoPrP89-230.pdb -d A 127 224 -r rlist-prp \
     -l ~/bin/eglib/fycqr.lib -o prp 1 1000 -maxc 20 -maxE 5000 -rlist -pdb
# library: /home/shbae/bin/eglib/fycqr.lib (21190)
# input PDB: MoPrP89-230.pdb (2209)
# rotation list: rlist-prp (43)
# testing steric clash of the input PDB coordinates
#       clash= 0 Evdw= 2779.36
```

The *rlist* file defines which residue to be randomly sampled for 
a set of phi, psi and chi dihedral angles, and which part
(N- or C-terminal from the residue) is to be rotated to 
make the selected phi, psi backbone dihedral angles.
For example, ```A 124 N``` means that phi,psi,chi dihedral angles are 
to be sampled and applied to the residue 124 by rotating residues 89-123.
Likewise, ```A 227 C``` means that dihedral angles are to be sampled 
and applied to the residue 227 by rotating residues 278-230.
As you might notice, all the residues before the rigid domain of interest should
be ```N``` and all the residues after the rigid domain of interest should 
be ```C```. This notation is necessary to define the *rigid domain of interest*
for proteins containing multiple rigid domains.

**eg** output shows the number of van der Waals clashes and 
van der Waals energy of the template pdb file. Any ensemble structure generated
by **eg** will have at least this number of clashes and this level of van der Waals energy
because coordinates of the rigid domain do not change in the ensemble structures.
You may adjust the maximum number of van der Waals clashes and 
maximum van der Waals energy by ```-maxc``` and ```-maxE```.
If these limits are set too low or the template ensemble structure is poorly defined, 
that is, junctions connecting rigid and disordered parts have constant clashes,
it will take long for **eg** to generate a desired number of ensemble structures 
that satisfy the criteria.

The above example would generate ```prp_0001.rot```,```prp_0001.pdb```,...,
```prp_1000.rot```,```prp_1000.pdb```. 

*.rot* files contain selected dihedral angles. 
In order to save disk space, you may delete the *.pdb* files and keep only the *.rot* files. 
You can rebuild the *.pdb* files using the *.rot* files and the initial template without any loss.

```
$ eg -i MoPrP89-230.pdb -r prp_0001.rot -rebuild
$ eg -i MoPrP89-230.pdb -r prp_0002.rot -rebuild
```

### Generate triangular surfaces

You need two sets of molecular surfaces or boundary elements: a static surface
for the rigid domain and an instantaneous surface for a snapshot of an ensemble structure.
The static surface serves as a reference boundary to which the instantaneous surface is
related to calculate the *velocity correlation*.
A matrix for the static surface does not undergo an inverse operation, 
so coarsening is not necessary.
In the mouse Prion(89-230), residues 127-224 can be assumed as rigid.
Thus a static surface is generated for the residues 127-224 and a series of 
instantaneous surfaces are generated for the whole molecule, the residues 89-230.

#### Static surface

```
$ pdb_to_xyzr MoPrP127-224.pdb | awk '{print $1,$2,$3,$4+1.1}' > MoPrP127-224.xyzr
$ msms.i86Linux2.2.6.1 -density 1 -if MoPrP127-224.xyzr -of MoPrP127-224 >& msms.log
$ msms2gts MoPrP127-224
```

One static surface is used with all instantaneous surfaces.

#### Instantaneous surface

```
$ pdb_to_xyzr prp_0001.pdb | awk '{print $1,$2,$3,$4+1.1}' > prp_0001.xyzr
$ msms.i86Linux2.2.6.1 -density 1 -if prp_0001.xyzr -of prp_0001 >& msms.log
$ msms2gts prp_0001 800
```

A coarsened instantaneous surface is genereated for each ensemble structure.

### Calculate diffusion tensors

```
$ be2 -gts prp_0001.0800 -rR MoPrP127-224 > prp_0001.out
$ cat prp_0001.out
Temperature  : 293.15 K
Viscosity    : 1.002 cP
Algorithm    : 1
GTS_surface  : prp_0001.0800.gts (Vertex: 402 Edge: 1200 Face: 800)
Dimension    : X   77.7835  Y   75.2813 Z   36.7675
GTS_surface  : MoPrP127-224.gts (Vertex: 6493 Edge: 19473 Face: 12982)
Dimension    : X   48.1480  Y   41.8170 Z   31.1910
Calculating G matrix (2400 x 2400) ...... time 00:00:24
Inverting   G matrix (2400 x 2400) ...... time 00:01:20
Surface_Area (  2.892479 ...  36.501948) Sum= 10729.6 A^2

Center_of_Diffusion 1.5827 6.6625 -0.4895
Dtt Eigenvalues (10^-7 cm^2 s^-1) | Eigenvectors
Dtt 1      9.646 |  6.712e-02 -1.479e-01  9.867e-01
Dtt 2     10.839 |  9.560e-01 -2.737e-01 -1.061e-01
Dtt 3     11.457 |  2.857e-01  9.504e-01  1.230e-01
Drr Eigenvalues (10^7 s^-1) | Eigenvectors
Drr 1      0.903 |  7.212e-02 -8.502e-02  9.938e-01
Drr 2      1.037 |  8.848e-01 -4.545e-01 -1.031e-01
Drr 3      1.353 |  4.604e-01  8.867e-01  4.244e-02
# gamma 6 eps 22 1/(6Diso)  15.186 ns
```

If you run BE2 for 10 different structures in an ensemble,
they would result in a range of rotational correlation times (1/(6Diso)).

```
$ grep Diso prp_*.out
prp_0001.out:# gamma 6 eps 22 1/(6Diso)  15.186 ns
prp_0002.out:# gamma 6 eps 22 1/(6Diso)  12.504 ns
prp_0003.out:# gamma 6 eps 22 1/(6Diso)  12.470 ns
prp_0004.out:# gamma 6 eps 22 1/(6Diso)  12.580 ns
prp_0005.out:# gamma 6 eps 22 1/(6Diso)  13.991 ns
prp_0006.out:# gamma 6 eps 22 1/(6Diso)  12.952 ns
prp_0007.out:# gamma 6 eps 22 1/(6Diso)  12.424 ns
prp_0008.out:# gamma 6 eps 22 1/(6Diso)  11.894 ns
prp_0009.out:# gamma 6 eps 22 1/(6Diso)  13.120 ns
prp_0010.out:# gamma 6 eps 22 1/(6Diso)  12.661 ns
```

For the mouse Prion(89-230), the experimental rotational correlation
time (1/(6Diso)) or tauc is 13.4 ns at standard condition (20 deg. C). 
Even for 10 structures, averaged tauc is pretty close to the experimental value.
However, a large number of ensemble structures are required for more
reproducible prediction. In the mouse Prion(89-230), it was found that 
averaging of about 1000-2000 ensemble structures leads to 
practical convergency at which tauc fluctuation is less than 0.2 ns.

### Temperature and viscosity

Hydrodynamic properties such as the translational and rotational 
diffusion tensors depend on the solution temperature and viscosity. 
These parameters can be specified by ```-t``` and ```-v``` options in **BE2**. 
By default, temperature is 293.13 K or 20 deg. C and viscosity is 1.002 cP.
A tauc value measured at different temperature can be converted to
a value at 20 deg. C by assuming that the buffer follows 
temperature dependent viscosity of water.

**Table for temperature dependent water viscosity**
```
#!/usr/bin/perl
print "# viscosity of water (cP)\n";
print "# (K) (deg) (cP)\n";
for ($C=15; $C < 40; $C++) 
{
    $K = $C + 273.15;
    printf("%.3f %.3f %.3f\n",$K,$C,2.414e-2*exp(247.8/($K-140.0)*log(10)));
}
```

**Convert tauc to the standard condition (water, 20 deg.C)**
```
#!/usr/bin/perl
# ex. conversion of tauc=8.4 at 25 deg.C to the standard condition
$tc = 8.4;
$C = 25.0;
$K = $C + 273.15;
$v = 2.414e-2*exp(247.8/($K-140.0)*log(10));
$v_= 2.414e-2*exp(247.8/(293.15-140.0)*log(10));
printf("tc: %f at %f C %f K --> tc: %f at 20 C\n",$tc,$C,$K,$tc*$K/293.15*$v_/$v);
```

### Comparison to the Stokes-Einstein estimation and false rigid assumption

For comparison, tauc can be estimated from the Stokes-Einstein equation
or tauc can be predicted following the above rigid BE2 method
with a false assumption of each structure in the ensemble being rigid.

**Stokes-Einstein estimation**

```
#!/usr/bin/perl
# ex. estimate tauc for mw=12719
$MW = 12719;
$SE = 1.38*6.022*293.15/1.002/6/($MW*0.73);
printf("MW:%d SE: %.5f (1/us) SE_predicted_tc: %.2f ns\n",$MW,$SE,1/(3*$SE));
```
