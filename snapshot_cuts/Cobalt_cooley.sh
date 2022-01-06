#!/bin/sh
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 12))

for i in 499 
#487 421 382 315 272 247 213 194 171 163
do 
echo $i

#mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/simulations_misc/snapshot_cuts/downsample /eagle/LastJourney/heitmann/OuterRim/M000/L4225/HACC000/analysis/Particles/STEP${i}/m000.mpicosmo.${i} /eagle/LastJourney/prlarsen/dc2_tests/bias/parts_${i}_reduced 100

mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/simulations_misc/snapshot_cuts/downsample /eagle/LastJourney/heitmann/ChristmasRuns/Ferenginar/HACC001/analysis/Particles/STEP499/m002.mpicosmo.499 /eagle/LastJourney/prlarsen/lj_bias/fe_parts_499_reduced_100 100

mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/simulations_misc/snapshot_cuts/downsample /eagle/LastJourney/heitmann/ChristmasRuns/Vulcan/HACC001/analysis/Particles/STEP499/m001.mpicosmo.499 /eagle/LastJourney/prlarsen/lj_bias/vu_parts_499_reduced_100 100

#/eagle/LastJourney/heitmann/LastJourney/M000P/L5025/HACC000/analysis/Particles/STEP${i}/m000p.mpicosmo.${i} /eagle/LastJourney/prlarsen/lj_bias/parts_${i}_reduced 100


#mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/simulations_misc/snapshot_cuts/downsample /eagle/LastJourney/heitmann/LastJourney/M000P/L5025/HACC000/analysis/Particles/STEP${i}/m000p.mpicosmo.${i} /eagle/LastJourney/prlarsen/lj_bias/parts_${i}_reduced 100

#mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/simulations_misc/snapshot_cuts/downsample_halo /eagle/LastJourney/heitmann/OuterRim/M000/L4225/HACC000/analysis/Halos/HaloCatalog/03_31_2018.OR.${i}.fofproperties /eagle/LastJourney/prlarsen/dc2_tests/bias/halos_${i}_reduced 1 1.e12 1.e15 
done
