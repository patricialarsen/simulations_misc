CC = ${HACC_MPI_CC}
CXX = ${HACC_MPI_CXX}

HPXDIR = /home/nfrontie/maps/Healpix_3.82
CFITSDIR = /home/nfrontie/maps/cfitsio-4.4.0/build

GIODIR = /home/nfrontie/trunk/HACC/aurora.kokkos_sycl/mpi
HACCDIR = /home/nfrontie/trunk/HACC


SCDIR = src

HPXCINCDIR = -I${HPXDIR}/include
HPXINCDIR = -I${HPXDIR}/include/healpix_cxx
HPXLIBDIR = ${HPXDIR}/lib
CFTLIBDIR = ${CFITSDIR}/lib 

GIOINCDIR = -I${HACCDIR}/submodules/genericio
GIOLIBDIR = ${GIODIR}/lib

HACCINCDIR = -I${HACCDIR}/cosmotools/algorithms/halofinder
HACCINCDIR2 = -I${HACCDIR}/cosmotools/common
HYDROINC = ${HACCDIR}/aurora.kokkos_sycl/mpi/nbody/common
HYDROINC2 = -I${HACCDIR}/nbody/common
HYDROINC3 = ${HACCDIR}/aurora.kokkos_sycl/mpi/cosmotools/algorithms/halofinder
BLOSCDIR = -I${GIODIR}/genericio/thirdparty/blosc


INCDIRS = ${GIOINCDIR} ${HPXINCDIR}  ${HPXCINCDIR}  ${BLOSCDIR} ${HACCINCDIR} ${HACCINCDIR2} ${HYDROINC2} -I${HYDROINC3} -I./inc/
LIBDIRS = -L${HPXLIBDIR}  -L${GIOLIBDIR} -L${CFTLIBDIR}

go_compile: ${SCDIR}/gravonly.cxx
	${CXX} -g -O3 -fopenmp -c  ${INCDIRS} ${SCDIR}/gravonly.cxx ${SCDIR}/utils_gravonly.cxx ${SCDIR}/pix_funcs_gravonly.cxx ${SCDIR}/PLParticles.cxx -std=c++17 ${HACC_MPI_CXXFLAGS}

go_link: ${SCDIR}/gravonly.cxx
	${CXX} -g  -O3 -fopenmp -o dens gravonly.o utils_gravonly.o pix_funcs_gravonly.o PLParticles.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++17 ${HACC_MPI_LDFLAGS}  


main_compile: ${SCDIR}/main.cxx 
	${CXX} -g -O3 -fopenmp -c  ${INCDIRS} ${SCDIR}/main.cxx ${SCDIR}/utils_dm.cxx ${SCDIR}/pix_funcs_dm.cxx ${SCDIR}/PLParticles.cxx -std=c++17 ${HACC_MPI_CXXFLAGS}

main_link: ${SCDIR}/main.cxx
	${CXX} -g  -O3 -fopenmp -o dens main.o utils_dm.o pix_funcs_dm.o PLParticles.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++17 ${HACC_MPI_LDFLAGS}  
	

hydro_compile: ${SCDIR}/hydro.cxx
	${CXX} -g -O3 -fopenmp -c  ${INCDIRS} ${SCDIR}/hydro.cxx ${SCDIR}/utils.cxx ${SCDIR}/pix_funcs.cxx ${SCDIR}/PLParticles.cxx ${SCDIR}/PLHalos.cxx -std=c++17

hydro_link: ${SCDIR}/hydro.cxx
	${CXX} -g -O3 -fopenmp  -o hydro_xray hydro.o utils.o pix_funcs.o PLParticles.o PLHalos.o ${HYDROINC3}/RadiativeCooling.o ${HYDROINC}/Domain.o ${GIODIR}/cosmotools/algorithms/halofinder/Partition.o ${GIODIR}/cosmotools/algorithms/halofinder/dims.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI  -L${CFTLIBDIR} -lcfitsio  -L${HPXLIBDIR} -lchealpix -lhealpix_cxx  -std=c++17


