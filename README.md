# Understanding BlasterMaster  
[Written by Andreas Luttens, Carlsson lab, Uppsala University - 2019]  
   
This is an attempt at documenting the protein grid preparation, blastermaster, in the the DOCK 3.7.X package. This should accelerate debugging and the understanding of future researchers working with DOCK. The crystal structure I worked with is 4EIY.  
  
By running blastermaster.py with the verbose toggled on, I tried to figure out what was going on under the hood of the machine.  
   
**Caveat emptor!** I document the protocol that uses the addNohydrogens flag. The hydrogenated protein, rec.crg.pdb, is prepared by another step, protonate.py from $d37scripts. 
## Identifying bindingsite residues: blasterBindingSiteResidues
The first step in blastermaster is to identify bindingsite residues. This is based on both the protein and the crystallised ligand and uses distance and sequence information. The first step uses rec.pdb and xtal-lig.pdb. 
Command:
```bash
$DOCKBASE/proteins/filt/bin/filt < $DOCKBASE/proteins/defaults/filt.params > filter.log
```
filt is an interactive script that asks a bunch of questions. To answer this unsupervised, there is a filt.params in the defaults directory that is used to parse the answers to the questions filt asks. The results from the calculation are stored in rec.site. This file contains every residue that is in the bindingsite. The residue type, its residue number and chain identifier are written down. If you want to visualize what residues filt deemed bindingsite residues, you can run:
```bash
for i in $(cat rec.site | awk '{print $2}' | sed 's|A||g') ; do (echo -e $i'+\c') ; done
```
This should give you a string of residue numbers. Copy it except for the last '+' character and open a pymol session with your rec.pdb and xtal-lig.pdb. Type in the commandline:
```bash
select res $selection
show sticks, sele
```
This should highlight the bindingsite residues.  
  
## Creating a molecular surface: blasterMolSurf
The next step in the protein grid generation is creating a molecular surface using ms (written by Michael Connolly). First the radii used by dms are copied from the defaults directory. The radii file contains atomic radii used by dms. The surface represents the Van der Waals surface of the protein. A grep command filters out lines from rec.pdb that contain the pattern HOH. dms crashes if it finds waters. It creates an intermediate file, rec.crg.pdb.dms, that serves as input for dms, command:
```bash
grep -a -v HOH rec.crg.pdb > rec.crg.pdb.dms
```
The same command is run on rec.site, to make sure no waters are marked as bindingsite residues.
```bash
grep -a -v HOH rec.site > rec.site.dms
```
Now the waterfree system is known, the following command is run if you run a default blastermaster calculation:
```bash
$DOCKBASE/proteins/dms/bin/dms rec.crg.pdb.dms -a -d 1.0 -i rec.site.dms -g dms.log -p -n -o rec.ms
```
As mentioned before, dms is the binary executable that calculates a molecular surface based on the protein and its bindingsite residues. The bindingsite file is used to tell dms to only generate the molecular surface around the bindingsite residues. Another parameter that determines what this surface will look like is the density parameter, the default is 1. The bigger the value, the denser the surface. It creates a logfile and a rec.ms, which is the molecular surface. The logfile mentions how many atoms were in your protein, how many surface points it created and what the density of points per area was eventually. To figure out what exactly this rec.ms file contains, you should look at the contents of the rec.site.dms file. Choose a protein residue number + chain identifier that is not in this list and grep for it in rec.ms:
```bash
grep " 8A" rec.ms
```
The lines that match are simply the coordinates of the atoms that are in this specific residue. Now grep for a residue that is in the rec.site.dms:
```bash
grep " 9A" rec.ms
```
Now more lines match the pattern. Not only the coordinates of the atoms in the bindingsite residue are written down, but coordinates of surface points too. If you examine the alpha-carbon of the residue TYR 9, CA, the first line has the atom coordinates, for instance: 
```bash
>TYR   9A   CA  -6.687    7.417    4.865 A
```
A couple of surface points for the same alpha-carbon are:
```bash
>TYR   9A   CA  -4.793    7.199    5.299 SR0  0.185  0.976 -0.163 -0.142
>TYR   9A   CA  -4.794    7.269    5.121 SS0  0.152  0.977 -0.213 -0.015
>TYR   9A   CA  -4.849    7.649    4.442 SC0  0.198  0.967  0.122 -0.223
```
For an atom record, the seventh field is "A." For a surface point record, the seventh field begins with an "S," followed by a "C," "R," or "S" according to whether the point is part of contact, reentrant, or saddle surface (saddle is a type of reentrant surface where the probe is in contact with exactly two atoms). This is followed by a digit used for depicting different density levels. The eighth field is the molecular surface area associated with the point in Å2. The last three fields are the unit normal vector pointing outward from the surface, which are written down because of the -n flag in ms. The normal unitvectors are used by sphgen in the next step.  
  
The rec.ms file can be opened with chimera, should you wish to visualize it.  
  
A rec.ts.ms is also created, which will be used for thin sphere selection. If you run a default blastermaster calculation, rec.ms and rec.ts.ms are identical.
## Creating initial spheres: blasterSphgen
The next step is to create an initial set of spheres. This calculation is based on the molecular surface built in the previous section and a set of parameters that are parsed in the INSPH file. The file contains the following lines:
```bash
>rec.ms
>R
>X
>0.
>5.0
>1.4
>all_spheres.sph
```
The parameters mean the following:
* rec.ms: molecular surface file built by ms
* R: sphere outside of surface (R) or inside surface (L)
* X: specifies subset of surface points to be used (X=all points)
* 0.0 prevents generation of large spheres with close surface contacts (default=0.0) 
* 5.0: maximum sphere radius in angstroms (default=4.0) 
* 1.4: minimum sphere radius in angstroms (default=radius of probe) 
* all_spheres.sph: clustered spheres file  
  
Sphgen generates sets of overlapping spheres to describe the shape of a molecule or molecular surface. For receptors, a negative image of the surface invaginations is created; for a ligand , the program creates a positive image of the entire molecule. Spheres are constructed using the molecular surface described by Richards (1977) calculated with the program dms . Each sphere touches the molecular surface at two points and has its radius along the surface normal of one of the points. For the receptor, each sphere center is outside the surface, and lies in the direction of a surface normal vector. For a ligand, each sphere center is inside the surface, and lies in the direction of a reversed surface normal vector.  
  
Spheres are calculated over the entire surface, producing approximately one sphere per surface point. This very dense representation is then filtered to keep only the largest sphere associated with each receptor surface atom. The filtered set is then clustered on the basis of radial overlap between the spheres using a single linkage algorithm. This creates a negative image of the receptor surface, where each invagination is characterized by a set of overlapping spheres. These sets, or clusters, are sorted according to numbers of constituent spheres, and written out in order of descending size. The largest cluster is typically the ligand binding site of the receptor molecule. The program showsphere writes out sphere center coordinates in PDB format and may be helpful for visualizing the clusters.  
  
The log OUTSPH contains a summary of the parsed parameters for sphgen.  
```bash
>density type = X
> reading  rec.ms                                                                             type   R
> # of atoms =   2797   # of surf pts =  12824
> finding spheres for   rec.ms                                                                          
> dotlim =     0.000
> radmax =    5.000
> Minimum radius of acceptable spheres?
>    1.400000    
> output to  all_spheres.sph                                                                 
> clustering is complete     10  clusters
 ```
 According to OUTSPH, there were 10 clusters found. A simple grep for the word clusters in all_spheres.sph, the outputfile, verifies this:
 ```bash
 grep clusters all_spheres.sph
 
>cluster     1   number of spheres in cluster    99
>cluster     2   number of spheres in cluster    15
>cluster     3   number of spheres in cluster     8
>cluster     4   number of spheres in cluster     6
>cluster     5   number of spheres in cluster     4
>cluster     6   number of spheres in cluster     4
>cluster     7   number of spheres in cluster     4
>cluster     8   number of spheres in cluster     2
>cluster     9   number of spheres in cluster     2
>cluster    10   number of spheres in cluster     2
>cluster     0   number of spheres in cluster  1296
 ```
The clusters are listed in numerical order from largest cluster found to the smallest. At the end of the clusters is cluster number 0. This is not an actual sphere cluster, but a list of all of the spheres generated whose radii were larger than the minimum radius, before the filtering heuristics (i.e., allowing only one sphere per atom and using a maximum radius cutoff) and clustering were performed.  
  
Looking into a few spheres of cluster 1, by running:
```bash
head all_spheres.sph  

>cluster     1   number of spheres in cluster    99
>   66  -1.64526   7.34550   3.26535   3.404 2409 0  0
>   69  -2.39692  10.29756   3.36131   3.179 2364 0  0
>   72  -6.66889  12.02877   8.40340   2.135  582 0  0
>  546  -6.90978   4.39463  17.55475   2.056  723 0  0
>  549  -5.74930   4.35000  16.89084   2.145  723 0  0
>  572  -5.54911  11.41246  15.70895   1.736 1496 0  0
>  575  -6.35453   4.81023  17.54406   2.247  723 0  0
>  580  -9.33744  13.31595   9.33673   1.653  583 0  0
>  581  -6.39474  17.25624  10.37895   3.177 2358 0  0
```
The values in the all_spheres.sph sphere file correspond to:  
  
* The number of the atom with which surface point i (used to generate the sphere) is associated.
* The x, y and z coordinates of the sphere center.
* The sphere radius.
* The number of the atom with which surface point j (second point used to generate the sphere) is associated.
* The critical cluster to which this sphere belongs.
* The sphere color. The color is simply an index into the color table that was specified in the header. Therefore, 1 corresponds to the first color in the header, 2 for the second, etc. 0 corresponds to unlabeled.  
  
**TODO: follow up**  
## Creating thin spheres: blasterThinSpheres
This step generates and selects thin spheres. The used scripts are thin_spheres.py and close_sph.py. The required input is the molecular surface that was generated earlier, rec.ts.ms. Several options can be parsed:
* ts_dist_ele: for low dielectric thin spheres, distance to protein surface, default = 1.0
* ts_radius_ele: for low dielectric thin spheres, radius of spheres, default = 1.0
* ts_dist_lds: for ligand desolvation thin spheres, distance to protein surface, default = 1.0
* ts_radius_lds: for ligand desolvation thin spheres, radius of spheres, default = 1.0
  
The outputfiles by default are called: low_die_thinspheres.sph & lig_des_thinspheres.sph. These files are identical if you run a default blastermaster calculation, since the default distances and radii are identical. The calculation command for the low dielectric thin spheres is the following: 
```bash
$DOCKBASE/proteins/thinspheres/thin_spheres.py -i rec.ts.ms -o low_die_thinspheres.sph -d 1.000000 -s 1.000000
```
**TODO: what if you want to visualize this step and create a pdb file**  
  
Now spheres with a parsable radius and distance from the protein surface were generated using a molecular surface that is based on the bindingsite residues. The following step is to select spheres that are close to the ligand, based on the xtal-lig.pdb file. The command which achieves this is the following:
```bash
python $DOCKBASE/proteins/thinspheres/close_sph.py low_die_thinspheres.sph xtal-lig.pdb low_die_thinspheres.sph.close 2.000000 1.000000
```
The inputfiles are the thinspheres, low_die_thinspheres.sph, from the previous calculation and the atom coordinates of the crystallised ligand, xtal-lig.pdb. The second to last argument is the threshold distance from the ligand. Only spheres closer than this threshold are stored in the outputfile, low_die_thinspheres.sph.close. The last argument is the sphere radius. You can verify there are less entries by the following command:
```bash
ls -ltr low_die_thinspheres.sph*
```
The .close file should contains less data than its parent inputfile. 
  
A similar (or identical, depending on the parsed radii) calculation is done for the ligand desolvation spheres.
## Converting ligand atoms into spheres: blasterPdbsph
The next step is to convert the heavy atoms of the crystallised ligand into spheres with a constant radius. This is achieved by running the following command:
```bash
$DOCKBASE/proteins/pdbtosph/bin/pdbtosph xtal-lig.pdb xtal-lig.match.sph
``` 
The radius of the heavy atoms is 0.7A. You can verify this by running:
```bash
cat xtal-lig.match.sph
```
The coordinates from the xtal-lig.pdb are preserved in this file.
**TODO: visualize the ligand spheres**
## Selecting low dielectric spheres: blasterLowDielectricSpheres
The next step is to select the low dielectric spheres. This calculation requires the heavy atom ligand spheres, xtal-lig.match.sph, the initial spheres that are based on the molecular surface of the bindingsite residues, all_spheres.sph, the hydrogenated receptor, rec.crg.pdb. The documentation calls rec.crg.pdb a charged receptor, this is because once the hydrogens have been assigned, protonation states are determined and therefore the protein structure will carry charges. A perl script will build low dielectric spheres based on the initial set of spheres and a minimum number of low dielectric spheres, hardcoded to 25, using the following command:
```bash
$DOCKBASE/proteins/makespheres1/makespheres1.cli.pl xtal-lig.match.sph all_spheres.sph rec.crg.pdb lowdielectric.sph 25 >& lowdielectric.spheres.log
``` 
A logfile, lowdielectric.spheres.log is created. The contents of this logfile look like this:
```bash
>There are 25 ligand heavy atoms
>Ligand center of coords (x y z): -0.42388 8.52728 17.12716
>useligsph flag does not exist.
>Starting number of spheres from all_spheres.sph cluster 0: 1296
>FINISHED READING IN FILE all_spheres.sph CLUSTER 0
>Number of sphere points after cutting off spheres too far from center of ligand (12 angstroms): 1054
>Number of spheres after removing too far (>7 angstroms) and too close to receptor (<1.2 angstroms): 1051
>Number of spheres polar: 269, nonpolar: 584, out of current total: 1051
>Extents of grid to search are X1 X2 Y1 Y2 Z1 Z2: -12.30531 11.49040 -3.43362 20.47639 5.18484 28.02714
>Number of spheres after removing spheres too close to each other (approximately <1.5 angstroms): 353
>Number of spheres after second pass removing spheres absolutely <0.8 from each other: 275
>No crystallographic spheres present, so determining sphere nearest center of ligand (x y z): -0.42388 8.52728 17.12716
>Center point: 96 xyz: 0.45782 8.60800 17.73916
>After continuity checking, number of spheres is: 215
>Spheres center of coords (x y z): -1.99754265116279 12.0964769302326 16.3856468372093
>Not crystallographic ligand, so keeping spheres closest to center of receptor, center of ligand, and center of spheres
>Receptor center of coords (x y z): -2.0203621737576 -8.58282767250626 17.7443496603504
>Average weighted-distance for final elimination is: 38.7597767582982
>Final number of output spheres is: 120
>Number of polar spheres: 39
>Number of nonpolar spheres: 67
>Number of spheres that are near neither polar nor nonpolar receptor atoms: 14
>Final number of spheres that are from crystallographic ligand: 0
``` 
Looking at the contents of the outputfile the perl script created:
```bash
head -n 20 lowdielectric.sph
>DOCK 5.2 ligand_atoms
>positive                       (1)
>negative                       (2)
>acceptor                       (3)
>donor                          (4)
>ester_o                        (5)
>amide_o                        (6)
>neutral                        (7)
>not_neutral                    (8)
>positive_or_donor              (9)
>negative_or_acceptor           (10)
>neutral_or_acceptor_or_donor   (11)
>donacc                         (12)
>cluster     1   number of spheres in cluster   120
> 9001   0.06400  16.91274  17.25932   3.044 1369 0  0
> 9002  -3.25654  16.23475  13.81292   2.400 2357 0  0
> 9003  -0.31266  16.71078  16.46934   3.462 1501 0  0
> 9004   1.08907  16.03707  14.54050   2.647 2333 0  0
> 9005 -10.29860  12.71042   9.59000   1.400  591 0  0
> 9006  -2.53163  16.29400  14.53131   2.622 2357 0  0
```
Then use grep for the first sphere entry coordinates on all_spheres.sph:
```bash
grep "0.06400  16.91274  17.25932"
>1501   0.06400  16.91274  17.25932   3.044 1369 0  0
```
The first field is an identifier, it can be used to convert the .sph file into a pdh and ensure unique residue indices. The next three fields are the coordinates of the sphere, which can be retrieved from the all_spheres.sph, the inputfile for this perl script. The last four entries are basically copied from the all_spheres.pdb and their meaning is described above. There is one field missing, and it is the initial atom i, which the sphere was based on in the inital sphere generation step described in the sections above.  
  
The spheres generated in this step will be used in several scripts downstream of the blastermaster procedure.  
## Converting spheres files into pdb's: blasterSphtopdb
During this step, three sets of .sph files are converted to .pdb's: the lowdielectric.sph, the low_die_thinspheres.sph.close and the lig_des_thinspheres.sph.close. This is done, for example, by the following command:
```bash
$DOCKBASE/proteins/showsphere/doshowsph.csh lowdielectric.sph 1 lowdielectric.sph.pdb >& lowdielectric.sph.pdb.log
```
Next to three .pdb files, three .log files will be created. The number 1 that is parsed to the cshell script denotes the cluster index identifier. The .pdb files can be opened with any structure viewer.
## Selecting matching spheres: blasterMatchingSpheres
In this step, the matching spheres will be determined by using the heavy atom ligand spheres, the initial set of spheres based on the molecular surface of the bindingsite residues and the hydrogenated protein structure. A perl script takes three extra parameters to create the matching_spheres.sph and a logfile using the following command:
```bash
$DOCKBASE/proteins/makespheres3/makespheres3.cli.pl 1.5 0.8 45 xtal-lig.match.sph all_spheres.sph rec.crg.pdb matching_spheres.sph >& matching_spheres.log
```
The three numerical parameters are:
* matching distance 1, corresponds with the gridsize, typical distance between two spheres: default = 1.5
* matching distance 2, corresponds with threshold that deems two spheres being too close: default = 0.8
* maximum number of matching spheres: default = 45  
Using grep to look for the word "cluster" verifies that there are at most 45 spheres.
```bash
grep cluster matching_spheres.sph
>cluster     1   number of spheres in cluster    45
```
The matching distance parameters are used to slice in the initial set of spheres that were derived from the molecular surface around bindingsite residues. The matching spheres are required to be proximate to the ligand, not to close to the protein surface, not too sparsely nor densely packed. You want a good representation of spheres that aligns with the dimensionality and discreteness of the grids we intend to build. Too little and you might miss out on valuable orientional sampling, too much and you will have to pay in computational resources for the redundant sampling.
  
The matching_spheres.sph can be eventually found in the dockfiles directory and will be used by DOCK for the orientational sampling.  
## Creating the box: blasterBox
The next step is to create the boundaries of the system of interest by constructing a box. A perl script achieves this by the following inputfiles, parameter and command:
```bash
$DOCKBASE/proteins/makebox/makebox.smallokay.pl xtal-lig.match.sph rec.crg.pdb box 10.0 >& makebox.log
```
The matching_spheres.sph are the heavy atom ligand spheres, rec.crg.pdb is the hydrogenated protein structure and the numerical parameter 10.0 stands for the box margin.  
  
You can check the content of the box file, it is in pdb format.
```bash
>HEADER    CORNERS OF BOX  -12.509  -8.586   4.541  11.566  26.456  29.777
>REMARK    CENTER (X Y Z)   -0.471   8.935  17.159
>REMARK    DIMENSIONS (X Y Z)   24.075  35.042  25.236
>ATOM      1  DUA BOX     1     -12.509  -8.586   4.541
>ATOM      2  DUB BOX     1      11.566  -8.586   4.541
>ATOM      3  DUC BOX     1      11.566  -8.586  29.777
>ATOM      4  DUD BOX     1     -12.509  -8.586  29.777
>ATOM      5  DUE BOX     1     -12.509  26.456   4.541
>ATOM      6  DUF BOX     1      11.566  26.456   4.541
>ATOM      7  DUG BOX     1      11.566  26.456  29.777
>ATOM      8  DUH BOX     1     -12.509  26.456  29.777
>CONECT    1    2    4    5
>CONECT    2    1    3    6
>CONECT    3    2    4    7
>CONECT    4    1    3    8
>CONECT    5    1    6    8
>CONECT    6    2    5    7
>CONECT    7    3    6    8
>CONECT    8    4    5    7
```
The dimension of the box will be used later on to trim the electrostatic grid.
## Concatenating inputfiles for grid construction: blasterCat
This step concatenates the hydrogenated protein structure and the low dielectric spheres to make an inputfile for Qnifft, a program that generated the electrostatic grids in a later step. The command that is used is:
```bash
cat rec.crg.pdb lowdielectric.sph.pdb > receptor.crg.lowdielectric.pdb
```
This is why the low dielectric spheres were converted into pdb format in a previous step.
## Creating the electrostatic grid: blasterQnifft
Blastermaster creates three types of grids: an electrostatic grid, a Van der Waals grid and a desolvation grid. Once the inputfiles concerning low dielectric spheres are available, Qnifft can construct the electrostatic grid. In the defaults directory there is a template for Qnifft calculations, called delphi.def. Qnifft requires an inputfile containing parameters. The entire template from the default is copied into a new file, named qnifft.parm. The following parameters are written into the qnifft.parm file:
* The gridsize: default = 193
* Amber charges: default = amb.crg.oxt
* Van der Waals radii: default = vdw.siz  
* pdb input: receptor.crg.lowdielectric.pdb  
  
The parameter file also contains the names of all the outputfiles:
* pdb output: qnifft.atm
* phi output: qnifft.electrostatics.phi

The following command constructs the electrostatic grid, together with logfiles.
```bash
$DOCKBASE/proteins/qnifft/bin/qnifft22_193_pgf_32 qnifft.parm >& qnifft.log
```
When we take a look at qnifft.atm:
```bash
head qnifft.atm

> GRASP PDB FILE
> FORMAT NUMBER= 1
> HEADER output from qnifft
> HEADER atom radii in columns 55-60
> HEADER atom charges in columns 61-67
>ATOM      1  C   ACE A   0      -5.110  20.966   1.421  1.90  0.5260    A    C
>ATOM      2  O   ACE A   0      -4.855  20.774   2.618  1.60 -0.5000    A    O
>ATOM      3  CH3 ACE A   0      -4.002  21.188   0.415  1.90 -0.0260    A    C
>ATOM      4  N   PRO A   1      -6.338  20.976   1.027  1.65 -0.2570    A    N
>ATOM      5  CA  PRO A   1      -7.497  20.769   1.930  1.86  0.1120    A    C
```
We see that the hydrogenated protein structure now also contains atom radii and partial amber charges.
  
After the electrostatic grid has been generated, we trim it using the dimension of the previously constructed box and the script, phiTrim.py:
```bash
python $DOCKBASE/proteins/blastermaster/phiTrim.py qnifft.electrostatics.phi box trim.electrostatics.phi
```
The reduced grid, trim.electrostatics.phi, can be eventually found in the dockfiles directory and will be used by DOCK to determine the strength of the electrostatic interactions using the partial charges of the atom types and their closest grid points.
## Creating the Van der Waals grid: blasterChemGrid
The next grid that will be generated is the Van der Waals grid. This is accomplished using the binary executable chemgrid. chemgrid reads an inputfile, INCHEM, containing parameters used in the grid construction. The contents of INCHEM look like:
```bash
>rec.crg.pdb
>prot.table.ambcrg.ambH
>vdw.parms.amb.mindock
>box
>0.2
>1
>4
>10
>2.3 2.6
>vdw
```
The parameters have the following meaning:
* rec.crg.pdb: hydrogenated protein structure
* prot.table.ambcrg.ambH: charge parameter file
* vdw.parms.amb.mindock: Van der Waals parameter file
* box: boundary conditions encapsulating system of interest
* 0.2: grid spacing in Ångstrom
* 1: electrostatic type; 1 distance-dependent dielectric , 2; constant dielectric
* 4: elecstrostatic scale for force field scoring
* 10: cutoff distance for energy calculations
* 2.3 2.8 ; bumping distances for polar and non-polar receptor atoms
* vdw: prefix of the output file -> vdw.vdw  
  
Chemgrid is used by the following command:
```bash
$DOCKBASE/proteins/chemgrid/bin/chemgrid >& vdw.log
```
A couple of outputfiles are generated, let's look at OUTCHEM first:
```bash
> receptor pdb file:
>rec.crg.pdb                                                                     
> receptor parameters will be read from:
>prot.table.ambcrg.ambH                                                          
> van der Waals parameter file:
>vdw.parms.amb.mindock                                                           
> input box file defining grid location:
>box                                                                             
> box center coordinates [x y z]:
>  -0.4710000       8.935000       17.15900    
> box x-dimension =    24.07500    
> box y-dimension =    35.04200    
> box z-dimension =    25.23600    
> grid spacing in angstroms
>   0.2000000    
> grid points per side [x y z]:
>          122         177         128
> total number of grid points =      2764032
> a distance-dependent dielectric will be used
>the dielectric function will be multiplied by   4.00
> cutoff distance for energy calculations:
>    10.00000    
> distances defining bumps with receptor atoms:
>receptor polar atoms  2.30
>receptor carbon atoms  2.60
> output grid prefix name:
>vdw  
``` 
Other files that were constructed are:
* vdw.bmp: bump grid 
* vdw.esp: electrostatic values for the receptor
* vdw.vdw : Van der Waals values for the receptor
* OUTPARM: messages pertaining to parameterization of receptor atoms; net charge on the receptor molecule including any ions or waters in the receptor pdb file
* PDBPARM: shows which parameters have been associated with each atom in the receptor pdb file  
  
There are three files related to Van der Waals interactions in the dockfiles directory. vdw.bmp is the bump grid, vdw.vdw are the Van der Waals values for the protein, vdw.parms.amb.mindock is used for rigid body minimization to optimize the best pose during live docking.
## Creating the desolvation grid: blasterSolvmap
Ligand desolvation is calculated as a sum of the atomic desolvation multiplied by a normalization factor that accounts for the extent to which the ligand atom is buried by the binding site. The atomic desolvation for each ligand atom can be calculated by AMSOL and is stored in the db2 file. The cost of desolvating each atom, or the normalization factor, is the distance weighted high dielectric volume displaced by the protein that is stored for each grid element in the active site.  
  
Two types of desolvation grids will be constructed, one for heavy atoms, and one for hydrogens. Two directories are therefore created, heavy and hydrogen. The following happens in both directories:  
  
First, the entire contents of the hydrogenated protein structure are copied into a file, rec.crg.lds.pdb, that will be used by solvmap. A temporary directory is created where a solvmap calculation will be submitted. Only one solvmap calculation per directory can be run, therefore, a temporary directory. A parameter file, INSEV, to control the solvmap calculation is constructed. The following is the content for the heavy/INSEV:
```bash
>rec.crg.lds.pdb
>ligand.desolv.heavy
>1.60,1.65,1.90,1.90,1.90,1.00
>1.4
>2
>box
>1.8
``` 
The parameters mean the following:
* the hydrogenated protein structure
* the desolvation map output
* radii of oxygen, nitrogen, carbon, sulfur, phosphorus, and "other" atoms.
* probe radius
* grid resolution
* box encapsulating system of interest
* assumed ligand atomic radius  
  
The following command starts the desolvation grid construction:
```bash
$DOCKBASE/proteins/solvmap/bin/solvmap >& solvmap.log
``` 
After solvmap is done constructing desolvation grids, both ligand.desolv.hydrogen and ligand.desolv.heavy are moved into the dockfiles directory. They will both be used to estimate the desolvation cost for placing an atom into a certain depth of the pocket. If said atom is a hydrogen, the hydrogen desolvation grid is used, otherwise DOCK uses the heavy desolvation grid.
## Combining blastermaster section, writing INDOCK: blasterIndock
The last step is the write the INDOCK file and string format a template file. For instance, if a default blastermaster calculation is run, the default INDOCK will be written, otherwise, names in the INDOCK file are altered. The phi gridsize is the only variable that is formatted and is inherited from the Qnifft step.
# How to contribute
This blastermaster.md serves as a small encyclopedia of what is roughly going on under the hood of blastermaster.py. It should lower the barrier for people to understand the procedure and help us steer away from black-box usage. Should new functionalities be introduced over time, or certain scripts change, then it is worthwhile to document this here. It is a small effort that can go a long way for the next generation in the lab.  
  
Cheers,  
Andreas
