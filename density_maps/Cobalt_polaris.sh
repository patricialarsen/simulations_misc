#!/bin/sh
#resoft
export OMP_NUM_THREADS=6
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/codes/healpix/Healpix_3.82/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/cfitsio-4.0.0/build/lib/
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 12))

for i in {500..624}
do 
  echo $i 
  mpirun -n 48 -f $COBALT_NODEFILE ./hydro /lus/eagle/projects/CosDiscover/nfrontiere/576MPC_RUNS/challenge_problem_576MPC_SEED_1.25e6_NPERH_AGN_2_NEWCOSMO/output/m000p.lc.mpicosmo.${i} /eagle/LastJourney/prlarsen/density_tests_2023/lc_dens_ne 2048 4 ${i} 0.7 1.0
done


