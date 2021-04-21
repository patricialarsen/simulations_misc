CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) -I/data/a/cpac/prlarsen/hacc/genericio -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/halofinder -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/mergertree_ct -I/data/a/cpac/prlarsen/hacc/cosmotools/common -I./src

#

OBJECTS := $(OBJDIR)/main.o $(OBJDIR)/Halos.o


#linking
match : $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/data/a/cpac/prlarsen/hacc/roomba.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o halo_check -fopenmp

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp
$(OBJDIR)/Halos.o: $(SRCDIR)/Halos.cxx
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/Halos.cxx -o $(OBJDIR)/Halos.o -fopenmp


clean:
	rm $(OBJDIR)/*.o halo_check
