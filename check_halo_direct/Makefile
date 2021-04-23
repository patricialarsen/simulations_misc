CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) -I/data/a/cpac/prlarsen/hacc/genericio -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/halofinder -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/mergertree_ct -I/data/a/cpac/prlarsen/hacc/cosmotools/common -I./src

#

OBJECTS := $(OBJDIR)/main.o $(OBJDIR)/Halos.o

OBJECTS2 := $(OBJDIR)/compare_dist.o $(OBJDIR)/Halos.o

#linking
match : $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/data/a/cpac/prlarsen/hacc/roomba.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o halo_check -fopenmp

dist : $(OBJECTS2) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS2) -L/data/a/cpac/prlarsen/hacc/roomba.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o compare_dist -fopenmp

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp
$(OBJDIR)/Halos.o: $(SRCDIR)/Halos.cxx
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/Halos.cxx -o $(OBJDIR)/Halos.o -fopenmp
$(OBJDIR)/compare_dist.o: $(SRCDIR)/compare_dist.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/compare_dist.cpp -o $(OBJDIR)/compare_dist.o -fopenmp


clean:
	rm $(OBJDIR)/*.o halo_check
