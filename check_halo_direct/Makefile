CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) -I/data/a/cpac/prlarsen/hacc/genericio -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/halofinder -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/mergertree_ct -I/data/a/cpac/prlarsen/hacc/cosmotools/common -I./src

#

OBJECTS := $(OBJDIR)/halo_check.o $(OBJDIR)/Halos_test.o

OBJECTS2 := $(OBJDIR)/compare_dist.o $(OBJDIR)/Halos_test.o

OBJECTS3 := $(OBJDIR)/match_pos.o $(OBJDIR)/Halos_test.o

#linking
halo_check: $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/data/a/cpac/prlarsen/hacc/roomba.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o halo_check -fopenmp

compare_dist : $(OBJECTS2) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS2) -L/data/a/cpac/prlarsen/hacc/roomba.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o compare_dist -fopenmp

match_pos : $(OBJECTS3) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS3) -L/data/a/cpac/prlarsen/hacc/roomba.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o match_pos -fopenmp


#compilation
$(OBJDIR)/halo_check.o: $(SRCDIR)/halo_check.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/halo_check.cpp -o $(OBJDIR)/halo_check.o -fopenmp
$(OBJDIR)/Halos_test.o: $(SRCDIR)/Halos_test.cxx
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/Halos_test.cxx -o $(OBJDIR)/Halos_test.o -fopenmp
$(OBJDIR)/compare_dist.o: $(SRCDIR)/compare_dist.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/compare_dist.cpp -o $(OBJDIR)/compare_dist.o -fopenmp
$(OBJDIR)/match_pos.o: $(SRCDIR)/nonid_match.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/nonid_match.cpp -o $(OBJDIR)/match_pos.o -fopenmp



clean:
	rm $(OBJDIR)/*.o halo_check compare_dist match_pos
