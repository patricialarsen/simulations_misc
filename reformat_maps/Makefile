CC = ${HACC_MPI_CC}
CXX = ${HACC_MPI_CXX}

ENVLOC = frontier.kokkos_hip

GIODIR = /ccs/home/prlarsen/codes/HACC/frontier.kokkos_hip/mpi
HACCDIR = /ccs/home/prlarsen/codes/HACC

HDF5FLAGS = -lhdf5
#-lhdf5_cpp -lhdf5
#-L/usr/lib64/  -lhdf5_cpp -lhdf5

SCDIR = src
HACCLIB = -lm -fopenmp

GIOINCDIR = -I${HACCDIR}/submodules/genericio
GIOLIBDIR = ${GIODIR}/lib

HDF5INCDIR = -I/sw/frontier/spack-envs/cpe23.12-cpu/opt/cce-17.0.0/hdf5-1.14.3-aaaha57nqzqzlc3bp2ldpwhhnyaxxts7/include
HDF5LIBDIR = -L/sw/frontier/spack-envs/cpe23.12-cpu/opt/cce-17.0.0/hdf5-1.14.3-aaaha57nqzqzlc3bp2ldpwhhnyaxxts7/lib

#HDF5INCDIR = -I/ccs/home/prlarsen/codes/HDF5/hdf5/build/include
#HDF5LIBDIR = -L/ccs/home/prlarsen/codes/HDF5/hdf5/build/lib
#HDF5INCDIR = -I/sw/frontier/spack-envs/cpe23.12-cpu/opt/cce-17.0.0/hdf5-1.14.3-aaaha57nqzqzlc3bp2ldpwhhnyaxxts7/include
#HDF5LIBDIR = -L/sw/frontier/spack-envs/cpe23.12-cpu/opt/cce-17.0.0/hdf5-1.14.3-aaaha57nqzqzlc3bp2ldpwhhnyaxxts7/lib

INCDIRS = ${GIOINCDIR}  -I./${SCDIR}/ ${HDF5INCDIR} -I./inc/
LIBDIRS =  -L${GIOLIBDIR} ${HDF5LIBDIR}


reformat_compile: ${SCDIR}/reformat_maps.cxx
	${CXX} -g -O3 -fopenmp -c  ${INCDIRS} ${SCDIR}/reformat_maps.cxx ${SCDIR}/utils_map.cxx ${SCDIR}/MapDef.cxx -std=c++17 ${HACC_MPI_CXXFLAGS}
reformat_link: ${SCDIR}/reformat_maps.cxx
	${CXX} -g  -O3 -fopenmp -o reformat_maps reformat_maps.o utils_map.o MapDef.o ${INCDIRS} -L${GIOLIBDIR}  -lGenericIOMPI -std=c++17 ${LIBDIRS} ${HACCLIB} ${HDF5FLAGS}



