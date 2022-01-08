# note: required export LD_LIBRARY_PATH=/home/prlarsen/usr/Healpix_3.31/lib/ at runtime

CC = mpicc
CXX = mpicxx

FLAGS = -I. 

HPXINCDIR = -I/home/prlarsen/usr/Healpix_3.80/include
HPXLIBDIR = -I/home/prlarsen/usr/Healpix_3.80/lib

HPXINCDIR2 = -I/home/prlarsen/usr/Healpix_3.80/include/healpix_cxx
FITSINCDIR = -I/home/prlarsen/usr/cfitsio-4.0.0/build/include
FITSLIBDIR = /home/prlarsen/usr/cfitsio-4.0.0/build/lib

GIOINCDIR = -I/home/prlarsen/hacc_fresh/trunk/genericio
BLOSCDIR = -I/home/prlarsen/hacc_fresh/trunk/genericio/thirdparty/blosc

GIOLIBDIR=/home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/lib
GIOINCDIR2= -I/home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio

#CXXHPXDIRS = -I/home/prlarsen/usr/Healpix_3.80/src/cxx/include
CXXHPXLIBS = -L/home/prlarsen/usr/Healpix_3.80/lib -lhealpix_cxx 
#-lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio


CFLAGS =  $(FITSINCDIR) $(HPXINCDIR) $(HPXINCDIR2) $(CXXHPXDIRS) $(GIOINCDIR)  $(GIOINCDIR2) $(BLOSCDIR) -fopenmp
LDFLAGS =  -lm $(CXXHPXLIBS)  -L$(HPXLIBDIR)  -L$(FITSLIBDIR)  -lGenericIOMPI  -L$(GIOLIBDIR) -lGenericIOMPI 
#-L/usr/lib/

SRCDIR = src
OBJDIR = obj
BINDIR = bin
#LIBDIR = usr/lib64/

cictogio: cic_gio2healpix.cxx
	$(CXX) -O3 -o cic_gio /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o ${HPXOFILES} ${CXXHPXDIRS} ${CFLAGS} ${LDFLAGS} ${FLAGS} cic_gio2healpix.cxx -lcurl 

cictocart: cic_gio2cart.cxx
	$(CXX) -O3 -o cic_cart /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o ${HPXOFILES} ${CXXHPXDIRS} ${CFLAGS} ${LDFLAGS} ${FLAGS} cic_gio2cart.cxx -lcurl 

#-dynamic

