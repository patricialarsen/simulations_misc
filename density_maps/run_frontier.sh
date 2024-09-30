#!/bin/bash
#SBATCH -A HEP142
#SBATCH -J lc_run
#SBATCH -o %x-%j.out
#SBATCH -t 00:30:00
#SBATCH -p batch
#SBATCH -N 256

source /ccs/home/prlarsen/codes/HACC/env/bashrc.frontier.cpu
cd /ccs/home/prlarsen/codes/simulations_misc/density_maps

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ccs/home/prlarsen/codes/Healpix_3.82/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ccs/home/prlarsen/codes/cfitsio-4.4.0/build/lib/
export OMP_NUM_THREADS=7

INPUT=/lustre/orion/hep114/proj-shared/Challenge_problems/Challenge_hydro/L1856_N7424/output/LC
#OUTPUT=/lustre/orion/hep142/proj-shared/prlarsen/maps/lc_tot
OUTPUT=/lustre/orion/hep142/proj-shared/Challenge_problems/Challenge_hydro/L1856_N7424/MAPS/
#lc_tot

# the 1 here is whether subfolders are turned on: 0=False, 1=True
# old
#srun -N 256 -n 2048 -c 7 ./hydro_xray ${INPUT}/ m000p.lc.mpicosmo. 1 ${OUTPUT} 8192 32 153 0.6766 1.0 153 5

srun -N 256 -n 2048 -c 7 ./hydro_xray ${INPUT}/ m000p.lc.mpicosmo. 1 ${OUTPUT} lc_tot 8192 32 153 0.6766 1.0 153 5 T 0.01

