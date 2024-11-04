CC = ${HACC_MPI_CC}
CXX = ${HACC_MPI_CXX}

#HPXDIR = /home/prlarsen/codes/healpix/Healpix_3.82_polaris
HPXDIR = /global/homes/p/plarsen/plarsen_git/Healpix_3.82
#CFITSDIR = /home/prlarsen/usr/cfitsio-4.0.0/build
CFITSDIR = /global/homes/p/plarsen/plarsen_git/cfitsio-4.4.0/build

GIODIR = /global/homes/p/plarsen/plarsen_git/HACC/perlmutter.kokkos_cuda/mpi
HACCDIR = /global/homes/p/plarsen/plarsen_git/HACC

SCDIR = src

#HPXDIR = /lustre/orion/hep114/scratch/prlarsen/Healpix_3.82
#CFITSDIR = /lustre/orion/hep114/scratch/prlarsen/cfitsio-4.2.0/build
#GIODIR = /lustre/orion/hep114/scratch/prlarsen/HACC/frontier.cpu/mpi
#HACCDIR = /lustre/orion/hep114/scratch/prlarsen/HACC

HPXCINCDIR = -I${HPXDIR}/include
HPXINCDIR = -I${HPXDIR}/include/healpix_cxx
HPXLIBDIR = ${HPXDIR}/lib
CFTLIBDIR = ${CFITSDIR}/lib 

GIOINCDIR = -I${HACCDIR}/submodules/genericio
GIOLIBDIR = ${GIODIR}/lib
HACCINCDIR = -I${HACCDIR}/cosmotools/algorithms/halofinder
HACCINCDIR2 = -I${HACCDIR}/cosmotools/common
HYDROINC = ${HACCDIR}/perlmutter.kokkos_cuda/mpi/nbody/common
HYDROINC2 = -I${HACCDIR}/nbody/common
HYDROINC3 = ${HACCDIR}/perlmutter.kokkos_cuda/mpi/cosmotools/algorithms/halofinder
BLOSCDIR = -I${GIODIR}/genericio/thirdparty/blosc


INCDIRS = ${GIOINCDIR} ${HPXINCDIR}  ${HPXCINCDIR}  ${BLOSCDIR} ${HACCINCDIR} ${HACCINCDIR2} ${HYDROINC2} -I${HYDROINC3} -I./inc/
LIBDIRS = -L${HPXLIBDIR}  -L${GIOLIBDIR} -L${CFTLIBDIR}



main_compile: main.cxx 
	${CXX} -c -fopenmp ${INCDIRS} ${SCDIR}/main.cxx utils.cxx pix_funcs.cxx PLParticles.cxx -std=c++11

main_link: main.cxx
	${CXX} -O3 -fopenmp -o dens main.o utils.o pix_funcs.o PLParticles.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++11
	

hydro_compile: ${SCDIR}/hydro.cxx
	${CXX} -c -fopenmp ${INCDIRS} ${SCDIR}/hydro.cxx ${SCDIR}/utils.cxx ${SCDIR}/pix_funcs.cxx ${SCDIR}/PLParticles.cxx ${SCDIR}/PLHalos.cxx -std=c++11

hydro_link: ${SCDIR}/hydro.cxx
	${CXX} -O3 -fopenmp -o hydro_xray hydro.o utils.o pix_funcs.o PLParticles.o PLHalos.o ${HYDROINC3}/RadiativeCooling.o ${HYDROINC}/Domain.o ${GIODIR}/cosmotools/algorithms/halofinder/Partition.o ${GIODIR}/cosmotools/algorithms/halofinder/dims.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++11


