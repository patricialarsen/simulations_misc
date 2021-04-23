## Halo direct comparison checks



The halo_check script directly compares two halo catalog outputs, assuming that these are one-to-one comparable. If the catalogs are not directly comparable (halo tags do not match or properties vary after evolution), you can instead use the compare_dist code for a first order check of the distributions.

### halo_check
The assumptions of the input files for halo_check are that the properties and halo tags match, however it does allow for different ordering between the files and for halo loss or addition. It also allows for name changes in the property files. 

After reading, re-ordering and correcting for missing halos, this computes the maximum differences in the quantities listed in halo_def_testing.h. If these ever exceed an input threshold, it computes the mean and standard deviation of the distribution of differences for each variable, as well as the maximum fractional difference and outputs these to the terminal. It then gives a summary of results with the number of variables failing the threshold, as well as the number of non-matching halos between the files. 

### compare_dist
Reads all the variables from each file, computing the mean and standard deviation. This then compares these values, and prints to terminal if these exceed a given threshold. It then gives a results summary with the number of variables failing the threshold as well as the total number of halos in the files. 



### How to use 

This should be built with the hacc environment, and linked to trunk. You should set OMP_NUM_THREADS to 1. 


In order to run this, alter the Makefile, type
```
make compare_dist
make halo_check
```
and run using
```

mpirun -n nprocesses ./compare_dist $FILE1 $FILE2 threshold
mpirun -n nprocesses ./halo_check $FILE1 $FILE2 threshold

```
a typical (truncated) output for halo_check is 

```

 sod_halo_angmom_x
 sod_halo_angmom_x
 ______________________________________
 mean difference = -1.08842e+11
 maximum fractional difference = 21436
 standard deviation of difference = 1.30895e+14
 mean of quantity = -4.05235e+10
 standard deviation of quantity = 2.19621e+14

 sod_halo_angmom_y
 sod_halo_angmom_y
 ______________________________________
 mean difference = -2.09428e+11
 maximum fractional difference = 9928.51
 standard deviation of difference = 1.35017e+14
 mean of quantity = -4.18144e+11
 standard deviation of quantity = 2.25129e+14

 sod_halo_angmom_z
 sod_halo_angmom_z
 ______________________________________
 mean difference = 8.32489e+09
 maximum fractional difference = 401576
 standard deviation of difference = 1.13304e+14
 mean of quantity = -1.39558e+11
 standard deviation of quantity = 1.73668e+14

 Results 
 ______________________________ 

 Comparison exceeded threshold of 0.01 for 7 variables
 out of a total of 51 variables 
 See above outputs for details  
 Total number of non-matching halos = 0



```
and a typical output for compare_dist is 
```

sod_halo_angmom_x
sod_halo_angmom_x
______________________________________
 means = -4.05235e+10 and = 6.8318e+10
 stddev = 2.19621e+14 and  = 1.62293e+14

sod_halo_angmom_y
sod_halo_angmom_y
______________________________________
 means = -4.18144e+11 and = -2.08716e+11
 stddev = 2.25129e+14 and  = 1.66669e+14

sod_halo_angmom_z
sod_halo_angmom_z
______________________________________
 means = -1.39558e+11 and = -1.47883e+11
 stddev = 1.73668e+14 and  = 1.34051e+14

 Results 
 ______________________________ 

 Comparison exceeded threshold of 0.01 for 6 variables
 out of a total of 51 variables 
 See above outputs for details  
 Difference in number of halos  = 0


```


### To-do
- Add integer options for comparisons

### Please note
This is largely based off Steve's match-up code to connect halo properties with halo lightcones, and re-uses code scripts from there. It also relies on the hacc code.
It is intended currently as a convenience script for internal hacc development. If you are using this to test hacc outupts and find that an iteration of the code has 
caused significant differences in the halo property files please contact me directly. 




## Python script 

This is a simple script reading in variables from both files and directly comparing outputs. Should only be used on small files.
