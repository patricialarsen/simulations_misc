CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) -I/home/prlarsen/hacc3/trunk/genericio -I/home/prlarsen/hacc3/trunk/cosmotools/algorithms/halofinder -I/home/prlarsen/hacc3/trunk/cosmotools/algorithms/mergertree_ct -I/home/prlarsen/hacc3/trunk/cosmotools/common
#

OBJECTS := $(OBJDIR)/main.o


#linking
make_hist : $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/home/prlarsen/hacc3/trunk/cooley.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o make_hist -fopenmp

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp

clean:
	rm $(OBJDIR)/*.o make_hist
