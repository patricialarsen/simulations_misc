CC = ${HACC_MPI_CC}
CXX = ${HACC_MPI_CXX}

ENVLOC = frontier.kokkos_hip

GIODIR = /ccs/home/prlarsen/codes/HACC/frontier.kokkos_hip/mpi
HACCDIR = /ccs/home/prlarsen/codes/HACC

SCDIR = src
HACCLIB = -lm -fopenmp

GIOINCDIR = -I${HACCDIR}/submodules/genericio
GIOLIBDIR = ${GIODIR}/lib


INCDIRS = ${GIOINCDIR}  -I./${SCDIR}/
LIBDIRS =  -L${GIOLIBDIR}

test_compile: ${SCDIR}/check_maps.cxx
	${CXX} -g -O3 -fopenmp -c  ${INCDIRS} ${SCDIR}/check_maps.cxx ${SCDIR}/map_def.cxx -std=c++17 ${HACC_MPI_CXXFLAGS}
test_link: ${SCDIR}/check_maps.cxx
	${CXX} -g  -O3 -fopenmp -o test_maps check_maps.o map_def.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI -std=c++17  ${HACCLIB}



