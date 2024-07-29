# note: required export LD_LIBRARY_PATH=/home/prlarsen/usr/Healpix_3.31/lib/ at runtime

CC = mpicc
CXX = mpicxx

FLAGS = -I. 

HACC_DIR = /home/prlarsen/hacc_fresh/HACC

HPXDIR = /home/prlarsen/codes/healpix/Healpix_3.82
HPXINCDIR = -I${HPXDIR}/include
HPXLIBDIR = ${HPXDIR}/lib
HPXINCDIR2 = -I${HPXDIR}/include/healpix_cxx

#HPXINCDIR = -I/home/prlarsen/usr/Healpix_3.80/include
#HPXLIBDIR = -I/home/prlarsen/usr/Healpix_3.80/lib

#HPXINCDIR2 = -I/home/prlarsen/usr/Healpix_3.80/include/healpix_cxx
FITSINCDIR = -I/home/prlarsen/usr/cfitsio-4.0.0/build/include
FITSLIBDIR = /home/prlarsen/usr/cfitsio-4.0.0/build/lib

GIOINCDIR = -I${HACC_DIR}/genericio
BLOSCDIR = -I${HACC_DIR}/genericio/thirdparty/blosc
BLOSCINCDIR = -I${HACC_DIR}/cooley.cpu/mpi/genericio/thirdparty/blosc

HACCINCDIR = -I${HACC_DIR}/cosmotools/algorithms/halofinder
HACCINCDIR2 = -I${HACC_DIR}/cosmotools/common

GIOLIBDIR=/home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/lib
GIOINCDIR2= -I/home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio


#CXXHPXLIBS = -L${HPXLIBDIR} -lhealpix_cxx
#CXXHPXDIRS = -I/home/prlarsen/usr/Healpix_3.80/src/cxx/include
#CXXHPXLIBS = -L/home/prlarsen/usr/Healpix_3.80/lib -lhealpix_cxx 
#-lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio


CFLAGS =  $(FITSINCDIR) $(HPXINCDIR) $(HPXINCDIR2) $(CXXHPXDIRS) $(GIOINCDIR)  $(GIOINCDIR2) $(BLOSCDIR) $(HACCINCDIR) $(HACCINCDIR2)  -fopenmp
LDFLAGS =  -lm $(CXXHPXLIBS)  -L$(HPXLIBDIR)  -L$(FITSLIBDIR)  -lGenericIOMPI  -L$(GIOLIBDIR) -lGenericIOMPI -lchealpix 
#-L/usr/lib/

SRCDIR = src
OBJDIR = obj
BINDIR = bin
#LIBDIR = usr/lib64/

HACC_INC_DIRS = ${GIOINCDIR} ${HACCINCDIR} ${HACCINCDIR2} ${BLOSCDIR} ${GIOINCDIR2}
HACC_LIB_DIRS = -L${GIOLIBDIR}

HPX_INC_DIRS = ${HPXINCDIR2} ${HPXINCDIR} ${FITSINCDIR}
HPX_LIB_DIRS = -L${HPXLIBDIR} -L${FITSLIBDIR}

cictogio: cic_gio2healpix.cxx
	$(CXX) -O3 -o cic_gio /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o ${HPXOFILES} ${CXXHPXDIRS} ${CFLAGS} ${LDFLAGS} ${FLAGS} cic_gio2healpix.cxx -lcurl 

cictocart: cic_gio2cart.cxx
	$(CXX) -O3 -o cic_cart /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o ${HPXOFILES} ${CXXHPXDIRS} ${CFLAGS} ${LDFLAGS} ${FLAGS} cic_gio2cart.cxx -lcurl 

memtest2: cic_gio2healpix_memoryupdate.cxx
	$(CXX) -O3 -o memtest /home/prlarsen/hacc_fresh/trunk/cooley.cpu/mpi/genericio/GenericIO.o ${HPXOFILES} ${CFLAGS} ${LDFLAGS} ${FLAGS} cic_gio2healpix_memoryupdate.cxx -lcurl 

memtest: cic_gio2healpix_memoryupdate.cxx
	$(CXX) -fopenmp -O3 -o memtest /home/prlarsen/hacc_fresh/HACC/cooley.cpu/mpi/genericio/GenericIO.o -lm -lGenericIOMPI ${HACC_INC_DIRS} ${HACC_LIB_DIRS} ${HPX_INC_DIRS} ${HPX_LIB_DIRS} -I${BLOSCINCDIR} -lchealpix cic_gio2healpix_memoryupdate.cxx -lchealpix -lcfitsio -lcurl 

memtest3: cic_gio2healpix_memoryupdate.cxx
	$(CXX) -fopenmp -O3 -o memtest /home/prlarsen/hacc_fresh/HACC/cooley.cpu/mpi/genericio/GenericIO.o -lm -lGenericIOMPI ${HACC_INC_DIRS} ${HACC_LIB_DIRS} ${HPX_INC_DIRS} ${HPX_LIB_DIRS} -I${BLOSCINCDIR} -lhealpix_cxx -lchealpix cic_gio2healpix_memoryupdate.cxx -lchealpix -lcfitsio -lcurl 


#-dynamic



# Add basic definitions 
#
