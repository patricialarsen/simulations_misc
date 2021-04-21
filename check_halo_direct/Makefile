CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) -I/home/prlarsen/hacc3/trunk/genericio -I/home/prlarsen/hacc3/trunk/cosmotools/algorithms/halofinder -I/home/prlarsen/hacc3/trunk/cosmotools/algorithms/mergertree_ct -I/home/prlarsen/hacc3/trunk/cosmotools/common -I./src
#

OBJECTS := $(OBJDIR)/main.o $(OBJDIR)/Halos.o


#linking
match : $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/home/prlarsen/hacc3/trunk/cooley.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o halo_check -fopenmp

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp
$(OBJDIR)/Halos.o: $(SRCDIR)/Halos.cxx
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/Halos.cxx -o $(OBJDIR)/Halos.o -fopenmp


clean:
	rm $(OBJDIR)/*.o halo_check
