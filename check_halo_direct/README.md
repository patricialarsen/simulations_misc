## Halo direct comparison checks


This testing suite directly compares two halo catalog outputs, assuming these are roughly one-to-one comparable. This inputs a parameters file params.txt, giving options for the tests performed

### Example params.txt script
```
/data/a/cpac/prlarsen/test_direc/arborx/final_may
/data/a/cpac/prlarsen/test_direc/arborx/finish
2
# 0 = ID matching, 1 = position matching, 2 = distribution comparisons, -1 = none
0
# 0 = SOD bin property file comparison (id matching), -1 = none
1
# 0 = bighaloparticle file matching (id matching), 1 = bhp and halopart and tags, -1 = none
256
# box size
0.0001
# threshold fraction
1.e13
# minimum mass for position matching
1.e15
# maximum mass for position matching


```
The first item here is the base analysis directory for the files to be read in. The second is the level of matching (ID-based, position-based or distribution only) for the halo catalogs. Any other number supresses the halo checks. The next option is whether or not to check the SOD bin property files, with ID based matching, and then the bighaloparticle files with ID-based matching, with an option to check the haloparticle and haloparticle tag files as well. This is currently just defaulting to the first subfile, as the memory requirements for a full check on these will be significant. 

The following inputs are 
- the box size of the simulation 
- a threshold for outputs (0 gives all outputs)
- a minimum mass for a position match
- a maximum mass for a position mass


After reading, re-ordering and correcting for missing halos, this computes the maximum differences in the quantities listed in halo_def_testing.h. If these ever exceed an input threshold, it computes the mean and standard deviation of the distribution of differences for each variable, as well as the maximum fractional difference and outputs these to the terminal. It then gives a summary of results with the number of variables failing the threshold, as well as the number of non-matching halos between the files. 


### How to use 

This should be built with the hacc environment, and linked to trunk. You should set OMP_NUM_THREADS to 1. 


In order to run this, alter the Makefile, type
```
make run_test
```
and run using
```

mpirun -n nprocesses ./run_test 

```
a typical (truncated) output is 

```

Partition 3D: [3:2:2]
 Parameter list
 /data/a/cpac/prlarsen/test_direc/arborx/final_may
 /data/a/cpac/prlarsen/test_direc/arborx/finish
 2
 0
 1
 256
 0.0001
 1e+13
 1e+15

 Comparing mean and standard deviations of property distributions
 ________________________________________________________________

Read 43 variables from /data/a/cpac/prlarsen/test_direc/arborx/final_may/m000p-499.haloproperties (6656625156 bytes) in 0.475155s: 13360.4 MB/s [excluding header read]
Read 43 variables from /data/a/cpac/prlarsen/test_direc/arborx/finish/m000p-499.haloproperties (6656625156 bytes) in 0.501408s: 12660.8 MB/s [excluding header read]

 Results 
 _______ 

 Comparison test passed! 
 All variables within threshold of 0.0001
 Difference in number of halos  = 0

 Checking consistency of sod halo bin files 
 ___________________________________________

Read 7 variables from /data/a/cpac/prlarsen/test_direc/arborx/final_may/m000p-delta200.0-499.sodpropertybins (1160980224 bytes) in 0.094238s: 11748.9 MB/s [excluding header read]
Read 7 variables from /data/a/cpac/prlarsen/test_direc/arborx/finish/m000p-delta200.0-499.sodpropertybins (1160980224 bytes) in 0.0930352s: 11900.8 MB/s [excluding header read]

 Results 
 _______ 

 Comparison test passed! 
 All variables within threshold of 0.0001
 Total number of non-matching halos = 0


 Checking constistency of halo tags and particle ids
 ___________________________________________________

Read 2 variables from /data/a/cpac/prlarsen/test_direc/arborx/final_may/m000p-499.bighaloparticles#0 (1419970864 bytes) in 1.10592s: 1224.49 MB/s [excluding header read]
Read 2 variables from /data/a/cpac/prlarsen/test_direc/arborx/finish/m000p-499.bighaloparticles#0 (1419970864 bytes) in 1.10761s: 1222.62 MB/s [excluding header read]

 Results 
 _______ 

 Total number of non-matching particles = 0


 Checking constistency of halo tags and particle ids
 ___________________________________________________

Read 2 variables from /data/a/cpac/prlarsen/test_direc/arborx/final_may/m000p-499.haloparticles#0 (2139138112 bytes) in 1.53598s: 1328.16 MB/s [excluding header read]
Read 2 variables from /data/a/cpac/prlarsen/test_direc/arborx/finish/m000p-499.haloparticles#0 (2139138112 bytes) in 1.48592s: 1372.91 MB/s [excluding header read]

 Results 
 _______ 

 Total number of non-matching particles = 0


 Checking constistency of halo tags and particle ids
 ___________________________________________________

Read 2 variables from /data/a/cpac/prlarsen/test_direc/arborx/final_may/m000p-499.haloparticletags#0 (2139138112 bytes) in 1.51808s: 1343.83 MB/s [excluding header read]
Read 2 variables from /data/a/cpac/prlarsen/test_direc/arborx/finish/m000p-499.haloparticletags#0 (2139138112 bytes) in 1.48677s: 1372.13 MB/s [excluding header read]

 Results 
 _______ 

 Total number of non-matching particles = 0

```


### To-do
- Add integer options for comparisons

### Please note
This is largely based off Steve's match-up code to connect halo properties with halo lightcones, and re-uses code scripts from there. It also relies on the hacc code.
It is intended currently as a convenience script for internal hacc development. If you are using this to test hacc outupts and find that an iteration of the code has 
caused significant differences in the halo property files please contact me directly. 




## Python script 

This is a simple script reading in variables from both files and directly comparing outputs. Should only be used on small files.
