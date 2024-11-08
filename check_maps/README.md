## Lightcone map checks 

This script runs basic validation tests on maps created on the lightcone particles. We assume for now an initial format of un-ordered GIO files with an ordering index, and both hydro and gravity-only outputs. The code is in parallel for more efficient reads. The checks that we run are:

- Confirm nside and number of pixels are as expected
- Confirm an output file exists for all steps 
- Output text file with the step-dependent mean, standard deviation, minimum, maximum and median of the maps
- Output text file with the number of zero-pixels per step
- Output text file with the number of NaN/inf pixels per step 
- Optionally create summed maps in five redshift slices for inspection to check map ordering and enable further validations


### How to use 

This should be built with the hacc environment, and linked to trunk.


In order to run this, alter the Makefile, type make and then run it with 
```
FOLDERS=T
HYDRO=T
PATH=/path/for/steps
MAP_NAME=/base_name_of_map
mpirun -n nprocesses ./map_check $PATH $FOLDERS $MAP_NAME $HYDRO
```


### Alterations
If the map code changes to output the maps with different names or extra quantities, update src/map_def.h

