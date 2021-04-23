CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) -I/data/a/cpac/prlarsen/hacc/genericio -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/halofinder -I/data/a/cpac/prlarsen/hacc/cosmotools/algorithms/mergertree_ct -I/data/a/cpac/prlarsen/hacc/cosmotools/common -I./src


OBJECTS := $(OBJDIR)/main.o $(OBJDIR)/LC_test.o


#linking
match : $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/data/a/cpac/prlarsen/hacc/roomba.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o lc_check -fopenmp

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp
$(OBJDIR)/LC_test.o: $(SRCDIR)/LC_test.cxx
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/LC_test.cxx -o $(OBJDIR)/LC_test.o -fopenmp

clean:
	rm $(OBJDIR)/*.o lc_check
