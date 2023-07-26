CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) -I/home/prlarsen/trunk_steve/trunk/genericio -I/home/prlarsen/trunk_steve/trunk/cosmotools/algorithms/halofinder -I/home/prlarsen/trunk_steve/trunk/cosmotools/algorithms/mergertree_ct -I/home/prlarsen/trunk_steve/trunk/cosmotools/common
#INCLUDES := -I $(SRCDIR) -I/home/prlarsen/hacc3/trunk/genericio -I/home/prlarsen/hacc3/trunk/cosmotools/algorithms/halofinder -I/home/prlarsen/hacc3/trunk/cosmotools/algorithms/mergertree_ct -I/home/prlarsen/hacc3/trunk/cosmotools/common
OBJECTS := $(OBJDIR)/main.o
#LIBS := /home/prlarsen/hacc3/trunk/cooley.cpu/mpi/lib/libGenericIOMPI.a
#LIBS := /home/prlarsen/trunk_steve/trunk/cooley.cpu/mpi/lib/libGenericIOMPI.a


#linking
match : $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/home/prlarsen/trunk_steve/trunk/cooley.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o match -fopenmp
#
#	$(CPP) $(CPPFLAGS) $(OBJECTS) -L/home/prlarsen/hacc3/trunk/cooley.cpu/mpi/lib -lCosmoTools -lGenericIOMPI -o match -fopenmp

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp 
#$(LIBS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp

clean:
	rm $(OBJDIR)/*.o match
