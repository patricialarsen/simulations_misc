#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#if defined(_OPENMP) 
#include <omp.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <stdint.h>

#include <iostream>
#include <fstream>

#include "GenericIO.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>

#include <vector>

#include <climits>
#include <algorithm>

using namespace gio;
using namespace std;


#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t

int main(int argc, char ** argv) {


   MPI_Init(&argc,&argv);
   int root_process = 0;
   int rank, numranks;
   MPI_Comm_size(MPI_COMM_WORLD, &numranks);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   string FileName(argv[1]);
   string NewFileName(argv[2]);
   float cx,cy,cz;
   float radius;
   //fraction = atoi(argv[3]);
   cx = atof(argv[3]);
   cy = atof(argv[4]);
   cz = atof(argv[5]);
   radius = atof(argv[6]);

   vector<POSVEL_T> xx, yy, zz, x_1, y_1, z_1;

   assert(sizeof(ID_T) == 8);

   size_t Np = 0;
   unsigned Method = GenericIO::FileIOMPI;

   double PhysOrigin[3], PhysScale[3];


 {
  GenericIO GIO(  MPI_COMM_WORLD, FileName ,Method);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  MPI_Barrier(MPI_COMM_WORLD);

  Np = GIO.readNumElems();
  printf("Np = %lu \n", Np);
  GIO.readPhysOrigin(PhysOrigin);
  GIO.readPhysScale(PhysScale);


  xx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
  yy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
  zz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));

  GIO.addVariable("x", xx, true);
  GIO.addVariable("y", yy, true);
  GIO.addVariable("z", zz, true);

  GIO.readData();
} 

  if(Np!=0){

   xx.resize(Np);
   yy.resize(Np);
   zz.resize(Np);

   for (int i=0; i<Np; i++){
   
   double dist;
   double dx = (double)(xx[i]-cx);
   double dy = (double)(yy[i]-cy);
   double dz = (double)(zz[i]-cz);
   dist = dx*dx + dy*dy + dz*dz;
   double rad_double;
   rad_double = (double)(radius*radius);
   if (dist<=rad_double){
       x_1.push_back(xx[i]);
       y_1.push_back(yy[i]);
       z_1.push_back(zz[i]);


    }
   }
}

{
    int CoordFlagsX = GenericIO::VarIsPhysCoordX;
    int CoordFlagsY = GenericIO::VarIsPhysCoordY;
    int CoordFlagsZ = GenericIO::VarIsPhysCoordZ;

    int Np2 = x_1.size();

    GenericIO NewGIO(MPI_COMM_WORLD, NewFileName); // change this to output name 

   x_1.resize(Np2 + NewGIO.requestedExtraSpace()/sizeof(POSVEL_T));
   y_1.resize(Np2 + NewGIO.requestedExtraSpace()/sizeof(POSVEL_T));
   z_1.resize(Np2 + NewGIO.requestedExtraSpace()/sizeof(POSVEL_T));

    NewGIO.setNumElems(Np2);

    for (int d = 0; d < 3; ++d) {
      NewGIO.setPhysOrigin(PhysOrigin[d], d);
      NewGIO.setPhysScale(PhysScale[d], d);
    }
 
    NewGIO.addVariable("x",x_1, CoordFlagsX | GenericIO::VarHasExtraSpace);
    NewGIO.addVariable("y",y_1, CoordFlagsY | GenericIO::VarHasExtraSpace);
    NewGIO.addVariable("z",z_1, CoordFlagsZ | GenericIO::VarHasExtraSpace);

    NewGIO.write();
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
   
}


