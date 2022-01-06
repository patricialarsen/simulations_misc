# note: required export LD_LIBRARY_PATH=/home/prlarsen/usr/Healpix_3.31/lib/ at runtime
CC = mpicc
CXX = mpicxx


GIOINCDIR = -I/home/prlarsen/hacc_fresh/trunk/genericio
BLOSCDIR = -I/home/prlarsen/hacc_fresh/trunk/genericio/thirdparty/blosc

GIOLIBDIR=/home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/lib
GIOINCDIR2= -I/home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio

CFLAGS = -O2  $(GIOINCDIR)  $(GIOINCDIR2) $(BLOSCDIR) -fopenmp
LDFLAGS =  -lm    -lGenericIOMPI  -L$(GIOLIBDIR) -lGenericIOMPI

SRCDIR = src
OBJDIR = obj
BINDIR = bin
LIBDIR = bin/lib


SOURCES  := downsample_snapshot.cxx downsample_halo_snapshot.cxx

downsample_halo: downsample_halo_snapshot.cxx
	${CXX} -O3 -o downsample_halo /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o -std=c++11  ${CFLAGS} ${FLAGS} downsample_halo_snapshot.cxx ${LDFLAGS} -lcurl ${HDF5FLAGS}


downsample: downsample_snapshot.cxx
	${CXX} -O3 -o downsample /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o -std=c++11  ${CFLAGS} ${FLAGS} downsample_snapshot.cxx ${LDFLAGS} -lcurl ${HDF5FLAGS}


cluster: cutout_cluster.cxx
	${CXX} -O3 -o cluster /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o -std=c++11  ${CFLAGS} ${FLAGS} cutout_cluster.cxx ${LDFLAGS} -lcurl ${HDF5FLAGS}

