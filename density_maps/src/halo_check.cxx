#include <mpi.h>

#include <stdint.h>


#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>
#include "stdio.h"
#include <fstream>
#include <sstream>
#include <algorithm>

#include <vector>
#include <unordered_map>

// healpix includes
#include <healpix_base.h>
#include "healpix_map.h"
#include "healpix_tables.h"
#include "math_utils.h"
#include "vec3.h"


//genericio includes
#include "GenericIO.h"
#include "BasicDefinition.h"


// local includes
#include "PLParticles.h"
#include "PLHalos.h"
#include "utils.h"
#include "pix_funcs.h"


#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


using namespace std;
using namespace gio;

//TODO: Nick: thread safety, general pixelization
int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int commRank, commRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

  double t1, t2, t3, t4, t5, t6;
  t1 = MPI_Wtime();

  // arguments that need to be added 
  bool borgcube = false;
  bool adiabatic = false;
  string cloudypath = "/pscratch/sd/p/plarsen/cloudy/";
  // string cloudypath = "/lus/eagle/projects/CosDiscover/nfrontiere/576MPC_RUNS/challenge_problem_576MPC_SEED_1.25e6_NPERH_AGN_2_NEWCOSMO/output/";

  #ifndef HYBRID_SG
     if (commRank==0)
       cout<< "hybrid_sg not defined" << endl;
     if (not adiabatic) {
       if (commRank==0)
         fprintf(stderr,"SUBGRID requested but subgrid definitions not enabled\n", argv[0]);
       exit(-1);
     }
  #else
     if (commRank==0)
       cout<< "hybrid_sg defined"<<endl;
     if (adiabatic){
       if (commRank==0)
         cout << "adiabatic option enabled, overwriting subgrid def";
     }
  #endif


  if(argc != 4) {
     if (commRank==0){
     fprintf(stderr,"USAGE: %s <inputfile> <hval> <samplerate>  \n", argv[0]);
     }
     exit(-1);
  }

  char filename[512];
  strcpy(filename,argv[1]);
  float hval = atof(argv[2]);
  float samplerate = atof(argv[3]);
 
  PLParticles P; 
  P.Allocate(); 

  MPI_Barrier(MPI_COMM_WORLD);

  read_particles(&P, filename);

  int status = check_xray_halo( &P, hval, borgcube, adiabatic,  samplerate, cloudypath);

  P.Deallocate();

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}


