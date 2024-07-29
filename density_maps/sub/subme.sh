#!/bin/bash
# UG Section 2.5, page UG-24 Job Submission Options
# Add another # at the beginning of the line to comment out a line
# NOTE: adding a switch to the command line will override values in this file.

# These options are MANDATORY at ALCF; Your qsub will fail if you don't provide them.
#PBS -A HLRedshift
#PBS -l walltime=02:00:00
#PBS -l select=64:system=polaris:ncpus=32:ngpus=4
#PBS -l filesystems=home:eagle
#PBS -l place=scatter
#ncpus = number of cores per node (probably should fix at 32)
#ngpus = number of gpus per node (probably should fix at 4)

# Highly recommended 
# The first 15 characters of the job name are displayed in the qstat output:
#PBS -N SZ_reduction

# If you need a queue other than the default (uncomment to use)
#PBS -q prod

# Controlling the output of your application
# UG Sec 3.3 page UG-40 Managing Output and Error Files
# By default, PBS spools your output on the compute node and then uses scp to move it the
# destination directory after the job finishes.  Since we have globally mounted file systems
# it is highly recommended that you use the -k option to write directly to the destination
# the doe stands for direct, output, error
#PBS -k doe

# If you want to merge stdout and stderr, use the -j option
# oe=merge stdout/stderr to stdout, eo=merge stderr/stdout to stderr, n=don't merge
#PBS -j oe


# The rest is an example of how an MPI job might be set up
source /home/prlarsen/hacc_fresh/HACC/env/bashrc.polaris.cpu
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
cd /home/prlarsen/codes/simulations_misc_new/simulations_misc/density_maps

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS=8           # Number of MPI ranks per node
NDEPTH=4           # Number of hardware threads per rank, spacing between MPI ranks on a node
NTHREADS=4         # Number of OMP threads per rank, given to OMP_NUM_THREADS

NTOTRANKS=$(( NNODES * NRANKS ))

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS}  THREADS_PER_RANK=${NTHREADS}"
#source bashrc.polaris.kokkos_cuda

#mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} --cpu-bind depth -env OMP_NUM_THREADS=${NTHREADS} ./set_affinity_gpu_polaris.sh ./hacc_tpm params/indat.params -n

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/codes/healpix/Healpix_3.82_polaris/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/cfitsio-4.0.0/build/lib/

for i in {333..474}
do
  echo $i 
  mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} --cpu-bind depth -env OMP_NUM_THREADS=${NTHREADS} ./hydro_p /lus/eagle/projects/CosDiscover/nfrontiere/576MPC_RUNS/challenge_problem_576MPC_SEED_1.25e6_NPERH_AGN_2_NEWCOSMO/output/m000p.lc.mpicosmo.${i} /eagle/LastJourney/prlarsen/density_tests_2023/lc_dens_ne 4096 16 ${i} 0.6766 1.0
done




