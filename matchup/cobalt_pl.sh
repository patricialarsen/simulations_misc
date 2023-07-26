#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES*12))
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
# probably want to divide this based on memory requirements - 98 is step 42-43, 1 is step 487-499
for i in {35..36}
# 36.. 80
#36..98}
#{1..74}
#8}
do
echo $i
i2=`expr $i + 1`
echo $i2
step=$(tail -n+$i steps_lj.txt | head -n1)
step2=$(tail -n+$i2 steps_lj.txt | head -n1)
echo $step
echo $step2
#mkdir /projects/LastJourney/prlarsen/halo_lightcone_LJ/output/lcHalos$step2
mkdir /projects/LastJourney/prlarsen/halo_lightcone_LJ/matched_output2/lcHalos$step2
LC_IN_FILE=/projects/LastJourney/prlarsen/halo_lightcone_LJ/output/lcHalos$step2/lc_intrp_halos.$step2
LC_OUT_FILE=/projects/LastJourney/prlarsen/halo_lightcone_LJ/matched_output2/lcHalos$step2/lc_intrp_halos_matched.$step2
#FOF_FILE=/projects/LastJourney/rangel/LastJourney/MergerTrees_updated_004/m000p-${step}.treenodes
FOF_FILE=/projects/LastJourney/heitmann/L5025/HACC000/analysis/Halos/b0168/Props/STEP${step}/m000p-${step}.haloproperties
echo $LC_IN_FILE
echo $LC_OUT_FILE
echo $FOF_FILE

mpirun --env LD_PRELOAD=$DARSHAN_PRELOAD -f $COBALT_NODEFILE  -n $PROCS /home/prlarsen/lc_codes/rangel/matchup2/match  $LC_IN_FILE $FOF_FILE $LC_OUT_FILE 

done

