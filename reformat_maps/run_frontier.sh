#!/bin/bash
#SBATCH -A HEP142
#SBATCH -J map_reformat
#SBATCH -o %x-%j.out
#SBATCH -t 02:00:00
##SBATCH -p batch
#SBATCH -q debug
#SBATCH -N 32
source /ccs/home/prlarsen/codes/HACC/env/bashrc.frontier.kokkos_hip
module load hdf5/1.14.3-mpi
#source /ccs/home/prlarsen/codes/HACC/env/bashrc.frontier.cpu
#export OMP_NUM_THREADS=7
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/frontier/spack-envs/cpe23.12-cpu/opt/cce-17.0.0/hdf5-1.14.3-aaaha57nqzqzlc3bp2ldpwhhnyaxxts7/lib
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/ccs/home/prlarsen/codes/HDF5/hdf5/build/lib

#srun -N 256 -n 2048 -c 7 ./hydro_xray ${INPUT}/ m000p.lc.mpicosmo. 1 ${OUTPUT} 8192 32 153 0.6766 1.0 153 5
#srun -n 64 ./test_maps /lustre/orion/hep142/proj-shared/INCITE/hydro/maps/old_analysis/healpix_maps/ lc_tot 1 1 16384 181 190 /lustre/orion/hep142/proj-shared/prlarsen/map_test/map_sum_181_190.gio
srun -N32 -n512 -c1 --cpu-bind=threads --threads-per-core=1 -m block:cyclic ./reformat_maps /lustre/orion/hep142/proj-shared/INCITE/hydro/maps/analysis/healpix_maps/  lc_tot 1 1 16384 240 624 /lustre/orion/hep142/proj-shared/prlarsen/maps_hdf5_hydro

