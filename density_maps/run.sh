#!/bin/bash
#SBATCH -A HEP142
#SBATCH -J calib_run
#SBATCH -o %x-%j.out
#SBATCH -t 02:00:00
#SBATCH -p batch
##SBATCH -q debug
#SBATCH -N 8

source bashrc.frontier.kokkos_hip
export OMP_NUM_THREADS=7

srun -N8 -n64 -c7 --gpus-per-node=8 --gpu-bind=closest ./hacc_tpm params/indat.params -n
