# note: required export LD_LIBRARY_PATH=/home/prlarsen/usr/Healpix_3.31/lib/ at runtime
CC = mpicc
CXX = mpicxx

FLAGS = -I. #-I/home/prlarsen/smoothing_files/src 

HPXINCDIR = -I/home/prlarsen/usr/Healpix_3.31/include
HPXLIBDIR = /home/prlarsen/usr/Healpix_3.31/lib

HPXINCDIR2 = -I/home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/include
HPXLIBDIR2 = /home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/lib

FITSINCDIR = -I/home/prlarsen/usr/cfitsio/lib/include
FITSLIBDIR = /home/prlarsen/usr/cfitsio/lib/lib

HACCDIR = /home/prlarsen/hacc_fresh/2022_hacc/HACC/
#HACCDIR = /home/prlarsen/hacc_fresh/HACC/

GIOINCDIR = -I${HACCDIR}/genericio
BLOSCDIR = -I${HACCDIR}/genericio/thirdparty/blosc

#GIOINCDIR= -I/projects/DarkUniverse_esp/prlarsen/hacc/genericio
#BLOSCDIR = -I/projects/DarkUniverse_esp/prlarsen/hacc/genericio/thirdparty/blosc

GIOLIBDIR=${HACCDIR}/cooley.cpu/mpi/lib
GIOINCDIR2= -I${HACCDIR}/cooley.cpu/mpi/genericio

AVDIR = -I/home/prlarsen/hacc_fresh/2022_hacc/HACC/cosmotools/algorithms/mergertree_ct
#AVDIR = -I/home/prlarsen/hacc_fresh/HACC/cosmotools/algorithms/mergertree_ct

HPXOFILES= /home/prlarsen/usr/Healpix_3.31/src/cxx/build.generic_gcc/Healpix_cxx/healpix_base.o

CXXHPXDIRS = -I/home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/include
CXXHPXLIBS = -L/home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/lib -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio

CFLAGS = -O2  $(FITSINCDIR) $(HPXINCDIR) $(CXXHPXDIRS) $(GIOINCDIR)  $(GIOINCDIR2) $(BLOSCDIR) -fopenmp
LDFLAGS =  -lm   $(CXXHPXLIBS)  -L$(HPXLIBDIR)  -lchealpix  -lcfitsio -L$(FITSLIBDIR)  -lGenericIOMPI  -L$(GIOLIBDIR) -lGenericIOMPI #-lhdf5 -lhdf5_cpp

HDF5FLAGS = -L/usr/lib64/  -lhdf5_cpp -lhdf5 #-I/home/prlarsen/codes/dtk/dtk


SRCDIR = src
OBJDIR = obj
BINDIR = bin
LIBDIR = bin/lib


pixelize_LJ: pixelize_halos_LJ.cxx
	${CXX} -O3 -o pixelize_LJ ${HACCDIR}/cooley.cpu/mpi/genericio/GenericIO.o   ${HPXOFILES} ${CXXHPXDIRS} ${CFLAGS} ${FLAGS} pixelize_halos_LJ.cxx utils_pixel_LJ.cxx Halos_LJ.cxx  ${LDFLAGS} -lcurl ${HDF5FLAGS} ${AVDIR}


