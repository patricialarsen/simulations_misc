# note: required export LD_LIBRARY_PATH=/home/prlarsen/usr/Healpix_3.31/lib/ at runtime
CC = mpicc
CXX = mpicxx

FLAGS = -I. -I/home/prlarsen/smoothing_files/src 

HPXINCDIR = -I/home/prlarsen/usr/Healpix_3.31/include
HPXLIBDIR = /home/prlarsen/usr/Healpix_3.31/lib

HPXINCDIR2 = -I/home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/include
HPXLIBDIR2 = /home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/lib

FITSINCDIR = -I/home/prlarsen/usr/cfitsio/lib/include
FITSLIBDIR = /home/prlarsen/usr/cfitsio/lib/lib

GIOINCDIR= -I/projects/DarkUniverse_esp/prlarsen/hacc/genericio
BLOSCDIR = -I/projects/DarkUniverse_esp/prlarsen/hacc/genericio/thirdparty/blosc

GIOLIBDIR=/projects/DarkUniverse_esp/prlarsen/hacc/cooley/mpi/lib
GIOINCDIR2= -I/projects/DarkUniverse_esp/prlarsen/hacc/cooley/mpi/genericio

HPXOFILES= /home/prlarsen/usr/Healpix_3.31/src/cxx/build.generic_gcc/Healpix_cxx/healpix_base.o

CXXHPXDIRS = -I/home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/include
CXXHPXLIBS = -L/home/prlarsen/usr/Healpix_3.31/src/cxx/generic_gcc/lib -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio

CFLAGS = -O2  $(FITSINCDIR) $(HPXINCDIR) $(CXXHPXDIRS) $(GIOINCDIR)  $(GIOINCDIR2) $(BLOSCDIR) -fopenmp
LDFLAGS =  -lm   $(CXXHPXLIBS)  -L$(HPXLIBDIR)  -lchealpix  -lcfitsio -L$(FITSLIBDIR)  -lGenericIOMPI  -L$(GIOLIBDIR) -lGenericIOMPI

HDF5FLAGS = -L/usr/lib64/  -lhdf5_cpp -lhdf5 -I/projects/DarkUniverse_esp/prlarsen/dtk 

SRCDIR = src
OBJDIR = obj
BINDIR = bin
LIBDIR = bin/lib


pixelize_sod: pixelize_halos_sod.cxx
	${CXX} -O3 -o pixelize_sod /projects/DarkUniverse_esp/prlarsen/hacc/cooley/mpi/genericio/GenericIO.o /projects/DarkUniverse_esp/prlarsen/dtk/hdf5_util.o  ${HPXOFILES} ${CXXHPXDIRS} ${CFLAGS} ${FLAGS} pixelize_halos_sod.cxx utils_pixel.cxx ${LDFLAGS} -lcurl ${HDF5FLAGS}


