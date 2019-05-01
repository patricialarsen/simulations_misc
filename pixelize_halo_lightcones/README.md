Divides halo GIO files into hdf5 pixel files in redshift ranges

Dependencies include:
Healpix C++ library
DTK toolkit (https://github.com/dkorytov/dtk - only requires hdf5 header files)
GIO library

To run a test version of this (on my setup) use:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/Healpix_3.31/lib/
mpirun -np 2 ./pixelize_sod /projects/DarkUniverse_esp/jphollowed/outerRim/lightcone_halos_octant_matchup_sod/lcHalos 8 z_0_01

Warning: there are many hard-coded oddities, such as missing steps, the step values of OuterRim, shift to lower octant, and correction for the great triangle.

This is intended only for the creation of DC2 halo files, general use is possible but may need alterations. 
