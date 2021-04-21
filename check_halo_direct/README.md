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

```

### To-do
Update variables

### Please note
This is largely based off Steve's match-up code to connect halo properties with halo lightcones, and re-uses code scripts from there. It also relies on the hacc code.
It is intended currently as a convenience script for internal hacc development. If you are using this to test hacc outupts and find that an iteration of the code has 
caused significant differences in the halo property files please contact me directly. 




## Python script 

This is a simple script reading in variables from both files and directly comparing outputs. Should only be used on small files.
