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

// local includes
#include "PLParticles.h"
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

  double t1, t2, t3, t4, t5;
  t1 = MPI_Wtime();

  if(argc != 10) {
     if (commRank==0){
     fprintf(stderr,"USAGE: %s <inputfile> <outputfile> <nside> <nside_low> <step> <hval> <samplerate> <start_step> <nsteps>  \n", argv[0]);
     }
     exit(-1);
  }

  char filename[512];
  strcpy(filename,argv[1]);
  const char *mpiioName = filename;

  int64_t nside = atoi(argv[3]);
  int64_t nside_low = atoi(argv[4]);
  string stepnumber = argv[5];
  string outfile = argv[2];

  float hval = atof(argv[6]);
  float samplerate = atof(argv[7]);
  int start_step = atoi(argv[8]);
  int nsteps = atoi(argv[9]);


  double rank1 = log2(nside);
  double rank2 = log2(nside_low);
  int rank_diff = (int) (rank1-rank2);
  printf("rank diff = %d \n", rank_diff);
 
  PLParticles P; 
  P.Allocate(); 

  if (commRank==0){ 
    printf("commRanks = %d \n", commRanks);
  }

  vector<int> send_count(commRanks,0);
  
  Healpix_Ordering_Scheme ring = RING;
  Healpix_Ordering_Scheme nest = NEST;
  T_Healpix_Base<int64_t> map_hires;
  T_Healpix_Base<int> map_lores;

  map_hires = Healpix_Base2 (nside, ring, SET_NSIDE); // we set to ring so interpolation is faster
  map_lores = Healpix_Base (nside_low, nest, SET_NSIDE);
  int64_t npix_lores = map_lores.Npix();
  int64_t npix_hires = map_hires.Npix();
  unordered_map<int64_t, int64_t> ring_to_idx;
  //raing_to_idx.reserve();


#ifdef _OPENMP
   if (commRank==0){
     printf("OpenMP is defined \n");
   }
#endif

  //TODO: alter to create the input pixelization and pass them in regardless

  vector<int> lores_pix;
  vector<int> pix_list_oct; // not needed for octant=0
  // need to create lists for octant=1,2 routines, currently works for full sky 
  get_pix_list_rank(0, commRank, commRanks,  npix_lores,  pix_list_oct, pix_list_oct, lores_pix);

  if (commRank==0){
  printf("Created get_pix_list_rank  \n");
  }

  vector<int64_t> start_idx;
  vector<int64_t> end_idx;
  vector<int64_t> pix_nums_start, pix_nums_end;

  int64_t map_size = lores_pix.size()*pow(4,rank_diff);
  ring_to_idx.reserve(map_size);


  int64_t count = 0; 
  vector<float> rho, phi, vel; // should generalize this so we're initializing it for an array of maps or a single map 
  for (int ii=0;ii< lores_pix.size() ;++ii){
  int pix_val = lores_pix[ii];
  initialize_pixel(pix_val, map_lores, map_hires, rho, phi, vel, count, start_idx,end_idx,pix_nums_start, pix_nums_end, rank_diff, &ring_to_idx);
  
  
  //initialize_pixel_hydro(pix_val, map_lores, map_hires, rho, phi, vel, count, start_idx,end_idx,pix_nums_start, pix_nums_end, rank_diff, &ring_to_idx);

  }

  
  MPI_Barrier(MPI_COMM_WORLD);
  t3 = MPI_Wtime();
  if (commRank==0){
    printf( "Elapsed time for initialization is %f\n", t3 - t1 );
  }

  for (int jj=0;jj<nsteps;jj++){

  MPI_Barrier(MPI_COMM_WORLD);
  t3 = MPI_Wtime();

  char step[10*sizeof(char)];
  char mpiioName[512];

  sprintf(step,"%d",start_step+jj);
  strcpy(mpiioName,mpiioName_base);
  strcat(mpiioName,step);



  read_and_redistribute(mpiioName, commRanks, &P, map_lores, map_hires, rank_diff);

  t2 = MPI_Wtime();
  if (commRank==0){
    printf( "Elapsed time for read and redistribute is %f\n", t2 - t3 );
  }



  int status = assign_dm_cic(rho, phi, vel, &P, map_hires, pix_nums_start, pix_nums_end,start_idx, ring_to_idx, hval, samplerate);
  //int status = assign_sz_ngp(ksz,tsz, &P, map_hires, pix_nums_start, pix_nums_end,start_idx, ring_to_idx);


  MPI_Barrier(MPI_COMM_WORLD);
  t4 = MPI_Wtime();
  if (commRank==0){
    printf( "Elapsed time for CIC is %f\n", t4 - t2 );
  }


  write_files(outfile, stepnumber,start_idx, end_idx, pix_nums_start, rho, phi, vel, npix_hires);
  MPI_Barrier(MPI_COMM_WORLD);
  t5 = MPI_Wtime();
  if (commRank==0){
    printf( "Elapsed time for write is %f\n", t5 - t4 );
    printf( "Elapsed time for step is %f\n", t5 - t3 );

  }
  }

  P.Deallocate();


  t2 = MPI_Wtime();
  if (commRank==0){
  printf( "Elapsed time is %f\n", t2 - t1 );
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}


