#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=12
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/Healpix_3.31/lib/

mpirun -np 12 ./pixelize_LJ /eagle/LastJourney/prlarsen/halo_lightcone_LJ/matched_output/lcHalos 32 z_2_3

