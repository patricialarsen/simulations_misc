## Halo direct comparison checks



This script directly compares two halo catalog outputs, assuming that these are one-to-one comparable. If the catalogs are not one-to-one comparable you should swap to the histogram checking method.

The assumptions of the input files are that the properties and halo tags match, however it does allow for different ordering between the files and for halo loss or addition. 

After reading, re-ordering and correcting for missing halos, this computes the maximum differences in the quantities:

## update

- fof_halo_center_x
- fof_halo_center_y
- fof_halo_center_z
- fof_halo_com_vx
- fof_halo_com_vy
- fof_halo_com_vz

and outputs them to the terminal. If these are non-zero, it computes the mean and standard deviation 
of the distribution of differences for each variable, as well as the maximum fractional difference.


### How to use 

This should be built with the hacc environment, and linked to trunk.


In order to run this, alter the Makefile, type make and then run it with 
```
mpirun -n nprocesses ./halo_check $FILE1 $FILE2
```
a typical output is, e.g. 

```
fof_halo_mass
______________________________________
 mean difference = 0
 maximum fractional difference = 0
 standard deviation of difference = 0
 mean of quantity = 2.96247e+11
 standard deviation of quantity = 3.95643e+12

fof_halo_ke
______________________________________
 mean difference = 6.56153e+06
 maximum fractional difference = 2.4446e-06
 standard deviation of difference = 3.39499e+10
 mean of quantity = 9.68491e+16
 standard deviation of quantity = 4.60448e+18

fof_halo_center_x
______________________________________
 mean difference = -1.64744e-08
 maximum fractional difference = 0.000502153
 standard deviation of difference = 2.29718e-05
 mean of quantity = 127.54
 standard deviation of quantity = 72.4255

fof_halo_center_y
______________________________________
 mean difference = 2.73343e-08
 maximum fractional difference = 0.000324737
 standard deviation of difference = 3.46306e-05
 mean of quantity = 125.18
 standard deviation of quantity = 72.6656

fof_halo_center_z
______________________________________
 mean difference = 6.55633e-08
 maximum fractional difference = 0.000861824
 standard deviation of difference = 8.76335e-05
 mean of quantity = 122.127
 standard deviation of quantity = 72.2195

fof_halo_angmom_x
______________________________________
 mean difference = 1.02118e+11
 maximum fractional difference = 217188
 standard deviation of difference = 8.24155e+14
 mean of quantity = 1.51496e+11
 standard deviation of quantity = 1.05413e+15


```

### To-do
Add integer comparisons

### Please note
This is largely based off Steve's match-up code to connect halo properties with halo lightcones, and re-uses code scripts from there. It also relies on the hacc code.
It is intended currently as a convenience script for internal hacc development. If you are using this to test hacc outupts and find that an iteration of the code has 
caused significant differences in the halo property files please contact me directly. 




## Python script 

This is a simple script reading in variables from both files and directly comparing outputs. Should only be used on small files.
