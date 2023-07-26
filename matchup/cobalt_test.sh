#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES*12))


mpirun -f $COBALT_NODEFILE  -n $PROCS /home/rangel/matchup/match  /projects/LastJourney/prlarsen/halo_lightcone_LJ/output/lcHalos76/lc_intrp_halos.76 /projects/LastJourney/heitmann/L5025/HACC000/analysis/Halos/b0168/Props/STEP77/m000p-77.haloproperties /projects/LastJourney/prlarsen/halo_lightcone_LJ/matched_output/a.out



#/projects/LastJourney/rangel/LastJourney/MergerTrees_updated_004/m000p-148.treenodes#6539 


#mpirun -f $COBALT_NODEFILE  -n $PROCS /home/prlarsen/lc_codes/rangel/matchup2/match  /projects/LastJourney/prlarsen/halo_lightcone_LJ/output/lcHalos76/lc_intrp_halos.76 /projects/LastJourney/rangel/LastJourney/MergerTrees_updated_004/m000p-77.treenodes /projects/LastJourney/prlarsen/halo_lightcone_LJ/matched_output/a.out
