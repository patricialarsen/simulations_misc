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
#include "utils_gravonly.h"
#include "pix_funcs_gravonly.h"

#ifdef HACC_NVRAM
#include "HaccIOUtil.h"
#endif


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

#ifndef HACC_NVRAM   
  if(argc != 14) {
     if (commRank==0){
     fprintf(stderr,"USAGE: %s <inputpath> <inputfile> <folders> <outpath> <outputfile> <nside> <nside_low> <step> <samplerate> <start_step> <nsteps> <output_downsampled> <downsampling_rate> \n", argv[0]);
     }
     exit(-1);
  }
#else
  if(argc != 15) {
     if (commRank==0){
     fprintf(stderr,"USAGE: %s <inputpath> <inputfile> <folders> <outpath> <outputfile> <nside> <nside_low> <step> <samplerate> <start_step> <nsteps> <output_downsampled> <downsampling_rate> <m_outInterBase> \n", argv[0]);
     }
     exit(-1);
  }
#endif



  char filepath[512];
  strcpy(filepath,argv[1]);
  const char *mpiioName_base = filepath;

  char filename[512];
  strcpy(filename,argv[2]);
  const char *mpiioName_file = filename;

  bool folders = atoi(argv[3]); // should be 0 if false

  int64_t nside = atoi(argv[6]);
  int64_t nside_low = atoi(argv[7]);
  string stepnumber = argv[8];
  //string outfile = argv[4];


  string outpath = argv[4];
  string outfile = argv[5];
  outfile = outpath + outfile;


  float samplerate = atof(argv[9]);
  int start_step = atoi(argv[10]);
  int nsteps = atoi(argv[11]);

  char output_downsampled = argv[12][0]; // T or F
  float downsampling_rate = atof(argv[13]);

#ifdef HACC_NVRAM
  string m_outInterBase = argv[14];
  int m_localRank, m_RPN, m_gio_nvram_file;
  assert(m_outInterBase != "");
  RankPerNode(m_RPN, m_localRank);
  m_gio_nvram_file = FileRankAssignment(commRanks/m_RPN,false,true);//File number for this rank on this node
  HACCNodeIOInfo::initialize(commRank, m_localRank, m_RPN, m_gio_nvram_file, m_outInterBase);//Initialize HACCNodeIO
#endif



  double rank1 = log2(nside);
  double rank2 = log2(nside_low);
  int rank_diff = (int) (rank1-rank2);
  if (commRank==0) 
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
  int64_t map_size = lores_pix.size()*pow(4,rank_diff);
  ring_to_idx.reserve(map_size);


  if (commRank==0){
  printf("Created get_pix_list_rank  \n");
  }

  vector<int64_t> start_idx;
  vector<int64_t> end_idx;
  vector<int64_t> pix_nums_start, pix_nums_end;


  int64_t count = 0; 
  vector<double> rho, phi, vel; // should generalize this so we're initializing it for an array of maps or a single map 

  for (int ii=0;ii< lores_pix.size() ;++ii){
  int pix_val = lores_pix[ii];
  // initialize all pixels on rank
  //
  initialize_pixel(pix_val, map_lores, map_hires, rho, phi, vel, count, start_idx,end_idx,pix_nums_start, pix_nums_end, rank_diff, &ring_to_idx);
  // make sure this is retaining this correctly
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

  if (folders){
  strcpy(mpiioName, mpiioName_base);
  strcat(mpiioName, "step_");
  strcat(mpiioName, step);
  strcat(mpiioName, "/");
  strcat(mpiioName, mpiioName_file);
  strcat(mpiioName, step);	  
  }
  else{
  strcpy(mpiioName, mpiioName_base);
  strcat(mpiioName, mpiioName_file);
  strcat(mpiioName, step);
  }
  
  for (int ii=0; ii<lores_pix.size(); ++ii){  
  clear_pixel(start_idx[ii], rank_diff,  rho, phi, vel);  
  }
    
  if (output_downsampled == 'T'){
    string downsampled_path_step = modifyPath(outpath, "downsampled_particles", string(step));
    string downsampled_file = downsampled_path_step + "/" + step + "_downsampled.gio";
    read_and_redistribute(mpiioName, commRanks, &P, map_lores, map_hires, rank_diff, true, downsampling_rate, downsampled_file);
  }
  else{
      read_and_redistribute(mpiioName, commRanks, &P, map_lores, map_hires, rank_diff);
  }

  //read_and_redistribute(mpiioName, commRanks, &P, map_lores, map_hires, rank_diff);

  MPI_Barrier(MPI_COMM_WORLD);
  t4 = MPI_Wtime();
  if (commRank==0){
    printf( "Redistribute time is %f\n", t4 - t3 );
  }

  int status = assign_dm_cic(rho, phi, vel, &P, map_hires, pix_nums_start, pix_nums_end,start_idx, ring_to_idx, samplerate);


  MPI_Barrier(MPI_COMM_WORLD);
  t5 = MPI_Wtime();
  if (commRank==0){
    printf("CIC time is %f\n", t5 - t4 );
  }


  write_files(outfile, step,start_idx, end_idx, pix_nums_start, rho, phi, vel, npix_hires);

  MPI_Barrier(MPI_COMM_WORLD);
  t6 = MPI_Wtime();
  if (commRank==0){
    printf( "Write time is %f\n", t6 - t5 );
  }

  t2 = MPI_Wtime();
  if (commRank==0){
  printf( "Elapsed time for this step is %f\n", t2 - t3 );
  }


  }

  P.Deallocate();


  t2 = MPI_Wtime();
  if (commRank==0){
  printf( "Elapsed time is %f\n", t2 - t1 );
  }

#ifdef HACC_NVRAM
 HACCIOShutdown(); //<-add this tricia
#endif
 
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}


