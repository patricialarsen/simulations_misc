## Lightcone checks


This script directly compares two lightcone particle outputs, assuming that these are one-to-one comparable.
For example this applies to comparisons between insitu lightcones from hacc restarts using the same parameter files or to post-processed lightcones. 

The assumptions of the input files are that the ID and replication values of the particles match, however it does allow for different ordering between 
the files and for particle loss or addition. 

After reading, re-ordering and correcting for missing particles, this computes the maximum differences in the quantities:

- x
- y
- z
- vx
- vy
- vz
- a
- phi

and outputs them to the terminal. If these are non-zero, it computes the mean and standard deviation 
of the distribution of differences for each variable, as well as the maximum fractional difference.


### How to use 

This should be built with the hacc environment, and linked to trunk.


In order to run this, alter the Makefile, type make and then run it with 
```
mpirun -n nprocesses ./lc_check $FILE1 $FILE2
```
a typical output is, e.g. 

```
Maximum dx = 0.000312805
Maximum dy = 0.000349522
Maximum dz = 0.000244141
Maximum dvx = 6.21346
Maximum dvy = 5.59067
Maximum dvz = 5.92481
Maximum da = 1.19209e-07
Maximum dphi = 3472.75
x position
______________________________________
 mean difference = 2.0434e-11
 maximum fractional difference = 0.00797565
 standard deviation of difference = 8.35351e-07

y position
______________________________________
 mean difference = -2.59427e-11
 maximum fractional difference = 0.163667
 standard deviation of difference = 8.93907e-07

z position
______________________________________
 mean difference = -1.39556e-10
 maximum fractional difference = 0.00185858
 standard deviation of difference = 8.64954e-07

x velocity
______________________________________
 mean difference = 6.60167e-08
 maximum fractional difference = 162.672
 standard deviation of difference = 0.00911104

y velocity
______________________________________
 mean difference = 1.49321e-07
 maximum fractional difference = 190.149
 standard deviation of difference = 0.00990435

z velocity
______________________________________
 mean difference = 1.0316e-07
 maximum fractional difference = 1436.48
 standard deviation of difference = 0.00981004

scale factor
______________________________________
 mean difference = -3.99949e-14
 maximum fractional difference = 1.25441e-07
 standard deviation of difference = 1.66848e-09

potential
______________________________________
 mean difference = 0.00478928
 maximum fractional difference = 1.34403
 standard deviation of difference = 5.294

```

### To-do
I intend to extend this by

- Including code to compare lightcones that do not have a one-to-one matching at the distribution level 
- Computing the angular power spectrum of the lightcone to test convergence on that property
- Allow for checks of halo lightcone properties 

### Please note
This is largely based off Steve's match-up code to connect halo properties with halo lightcones, and re-uses code scripts from there. It also relies on the hacc code.
It is intended currently as a convenience script for internal hacc development. If you are using this to test hacc outupts and find that an iteration of the code has 
caused significant differences in the lightcones please contact me directly. 

Also please do note that non-zero differences at the level of machine precision (and somewhat higher after evolution from a restart) are expected and not cause for alarm.
These often cause a handful of particles to belong to different lightcone outputs and small number differences do not signify true particle loss or creation in the simulation. 




## Python script 

This is a more in-depth but less scalable testing code. It runs as a serial python script, reading in the data, making density plots, making sure missing ids are on the boundaries, creating histograms of values and plotting an angular power spectrum comparison. Once you have confirmed that differences exist using the lc_check code you may want to use this to look at the output comparisons in more detail. 
