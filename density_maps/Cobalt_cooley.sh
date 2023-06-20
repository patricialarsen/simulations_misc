#!/bin/sh
#resoft
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/Healpix_3.31/lib/
export GENERICIO_RANK_PARTITIONS=12
export GENERICIO_PARTITIONS_USE_NAME=0
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 12))


#mpirun -f $COBALT_NODEFILE -n $PROCS  /home/prlarsen/codes/density_maps/density_maps_git/density_maps/ngp /projects/LastJourney/heitmann/L5025/HACC000/analysis/Lightc_full/STEP421/m000p.lc.mpicosmo.421 8192 2 421


#/home/prlarsen/codes/density_maps/density_maps_git/density_maps


#mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/density_maps/density_maps/ngp /projects/LastJourney/heitmann/L5025/HACC000/analysis/Lightc_full/STEP421/m000p.lc.mpicosmo.421 8192 2 421


#for i in {324..499}
#do
#  echo $i
#  mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/density_maps/density_maps_git/density_maps/ngp /projects/LastJourney/heitmann/L5025/HACC000/analysis/Lightc_full/STEP${i}/m000p.lc.mpicosmo.${i} 8192 4 ${i}
#done


for i in {210..219}
#{235..307}
do
  echo $i
  mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/density_maps/density_maps_git/density_maps/ngp /projects/LastJourney/heitmann/L5025/HACC000/analysis/Lightc_full/STEP${i}/m000p.lc.mpicosmo.${i} 8192 8 ${i}
done

for i in {123..185}
#{235..307}
do
  echo $i
  mpirun -f $COBALT_NODEFILE -n $PROCS /home/prlarsen/codes/density_maps/density_maps_git/density_maps/ngp /projects/LastJourney/heitmann/L5025/HACC000/analysis/Lightc_full/STEP${i}/m000p.lc.mpicosmo.${i} 8192 8 ${i}
done

