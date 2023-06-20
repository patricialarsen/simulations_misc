# density_maps
Create density maps from GIO lightcone outputs

Run make ngp_compile 
then make ngp_link

then run, e.g. 
mpirun -np 4 ./ngp  /data/a/cpac/prlarsen/STEP421/m000p.lc.mpicosmo.421 8192 2

mpirun -np 4 ./ngp /projects/LastJourney/heitmann/L5025/HACC000/analysis/Lightc_full/STEP499/m000p.lc.mpicosmo.499 1024 2 

remember to write to a new file



source /projects/DarkUniverse_esp/prlarsen/hacc/env/bashrc.cooley.default
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/Healpix_3.31/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/genericio/frontend

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/tiff-4.0.9/install/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/CGAL-4.13/install/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/mpfr-4.0.1/install/lib/
export LD_LIBRARY_PATH


export OMP_NUM_THREADS=3
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/codes/healpix/Healpix_3.82/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prlarsen/usr/cfitsio-4.0.0/build/lib/

mpirun -n 4 ./dens /lus/eagle/projects/CosDiscover/nfrontiere/576MPC_RUNS/challenge_problem_576MPC_GRAV_ONLY/output/m000p.lc.mpicosmo.624 /eagle/LastJourney/prlarsen/density_tests_2023/lc_dens 2048 4 624
