CC = ${HACC_MPI_CC}
CXX = ${HACC_MPI_CXX}

HPXDIR = /home/prlarsen/codes/healpix/Healpix_3.82_polaris
CFITSDIR = /home/prlarsen/usr/cfitsio-4.0.0/build
GIODIR = /home/prlarsen/hacc_fresh/HACC/polaris.kokkos_cuda/mpi
HACCDIR = /home/prlarsen/hacc_fresh/HACC


#HPXDIR = /lustre/orion/hep114/scratch/prlarsen/Healpix_3.82
#CFITSDIR = /lustre/orion/hep114/scratch/prlarsen/cfitsio-4.2.0/build
#GIODIR = /lustre/orion/hep114/scratch/prlarsen/HACC/frontier.cpu/mpi
#HACCDIR = /lustre/orion/hep114/scratch/prlarsen/HACC

HPXCINCDIR = -I${HPXDIR}/include
HPXINCDIR = -I${HPXDIR}/include/healpix_cxx
HPXLIBDIR = ${HPXDIR}/lib
CFTLIBDIR = ${CFITSDIR}/lib 

GIOINCDIR = -I${HACCDIR}/genericio
GIOLIBDIR = ${GIODIR}/lib
HACCINCDIR = -I${HACCDIR}/cosmotools/algorithms/halofinder
HACCINCDIR2 = -I${HACCDIR}/cosmotools/common
HYDROINC = ${HACCDIR}/polaris.kokkos_cuda/mpi/nbody/common
HYDROINC2 = -I${HACCDIR}/nbody/common

BLOSCDIR = -I${GIODIR}/genericio/thirdparty/blosc


INCDIRS = ${GIOINCDIR} ${HPXINCDIR}  ${HPXCINCDIR}  ${BLOSCDIR} ${HACCINCDIR} ${HACCINCDIR2} ${HYDROINC2}
LIBDIRS = -L${HPXLIBDIR}  -L${GIOLIBDIR} -L${CFTLIBDIR}



main_compile: main.cxx 
	${CXX} -c -fopenmp ${INCDIRS} main.cxx utils.cxx pix_funcs.cxx PLParticles.cxx -std=c++11

main_link: main.cxx
	${CXX} -O3 -fopenmp -o dens main.o utils.o pix_funcs.o PLParticles.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++11
	

hydro_compile: hydro.cxx
	${CXX} -c -fopenmp ${INCDIRS} hydro.cxx utils.cxx pix_funcs.cxx PLParticles.cxx -std=c++11

hydro_link: hydro.cxx
	${CXX} -O3 -fopenmp -o hydro_p hydro.o utils.o pix_funcs.o PLParticles.o ${HYDROINC}/RadiativeCooling.o ${HYDROINC}/Domain.o /home/prlarsen/hacc_fresh/HACC/polaris.kokkos_cuda/mpi/cosmotools/algorithms/halofinder/Partition.o /home/prlarsen/hacc_fresh/HACC/polaris.kokkos_cuda/mpi/cosmotools/algorithms/halofinder/dims.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++11

bhp_compile: bhp.cxx
	${CXX} -c -fopenmp ${INCDIRS} bhp.cxx utils.cxx pix_funcs.cxx PLParticles.cxx PLHalos.cxx -std=c++11

bhp_link: bhp.cxx
	${CXX} -O3 -fopenmp -o adiabatic_bhp bhp.o utils.o pix_funcs.o PLParticles.o PLHalos.o ${HYDROINC}/RadiativeCooling.o ${HYDROINC}/Domain.o /home/prlarsen/hacc_fresh/HACC/polaris.kokkos_cuda/mpi/cosmotools/algorithms/halofinder/Partition.o /home/prlarsen/hacc_fresh/HACC/polaris.kokkos_cuda/mpi/cosmotools/algorithms/halofinder/dims.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++11

