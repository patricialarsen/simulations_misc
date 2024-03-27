#!/bin/bash
#SBATCH -A m4075
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 0:30:00
#SBATCH --nodes=4
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=16

source /opt/cray/pe/cpe/23.12/restore_lmod_system_defaults.sh
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

#shifter --image=lsstsqre/centos:7-stack-lsst_distrib-w_2024_10 /bin/bash
#setup lsst_distrib

source /global/homes/p/plarsen/plarsen_git/HACC/env/bashrc.perlmutter.kokkos_cuda
cd /global/homes/p/plarsen/plarsen_git/simulations_misc/density_maps

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/homes/p/plarsen/plarsen_git/Healpix_3.82/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/homes/p/plarsen/plarsen_git/cfitsio-4.4.0/build/lib/
export OMP_NUM_THREADS=8

srun -n 32 --cpu-bind=cores -c16 ./hydro_xray /pscratch/sd/p/plarsen/HACC_LC/590/m000p.lc.mpicosmo.  /pscratch/sd/p/plarsen/HACC_LC/590/lc_hydro_tot 4096 2 590 0.6766 1.0 590 10

