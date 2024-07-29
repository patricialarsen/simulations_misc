#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=6

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/cfitsio-4.0.0/build/lib/
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/Healpix_3.82/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/codes/healpix/Healpix_3.82/lib
OMP_NUM_THREADS=1
export OMP_NUM_THREADS

export LD_LIBRARY_PATH


mpirun -f $COBALT_NODEFILE -n $PROCS ./memtest /lus/eagle/projects/CosDiscover/nfrontiere/576MPC_RUNS/challenge_problem_576MPC_ADIABATIC/output/m000p.lc.mpicosmo. /eagle/LastJourney/prlarsen/LC_hydro/2023/ksz_310_400.fits /eagle/LastJourney/prlarsen/LC_hydro/2023/tsz_310_400.fits 0.01 2 600 2048

#8192

#correct_massfunc.py
