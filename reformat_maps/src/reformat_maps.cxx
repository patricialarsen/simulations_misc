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

#include <vector>
#include <climits>
#include <algorithm>

//#include "H5Cpp.h"

#include "GenericIO.h"
#include "MapDef.h"

#include "utils_map.h"

using namespace gio;
using namespace std;
//using namespace H5;

#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


int main(int argc, char ** argv) {

    MPI_Init(&argc,&argv);
    int root_process = 0;
    //int rank, numranks;
    int commRank, commRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &commRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
 
 // NEW
  if(argc != 9) {
     if (commRank==0){
     fprintf(stderr,"USAGE: %s <inputpath> <inputfile> <folders> <hydro> <nside> <step_start> <step_end> <outfile> \n", argv[0]);
     }
     exit(-1);
  }

  char filepath[512];
  char filename[512];
  //char outfile[512];
  strcpy(filepath, argv[1]);
  strcpy(filename, argv[2]);
  //strcpy(outfile, argv[8]);
  string outfile = argv[8];

  const char *mpiioName_base = filepath;
  const char *mpiioName_file = filename;

  bool folders = atoi(argv[3]); // should be 0 if false
  bool hydro = atoi(argv[4]); // should be 0 if false
  int64_t nside = atoi(argv[5]);
  int step_start = atoi(argv[6]);
  int step_end = atoi(argv[7]);

  //
  int64_t nside_long = (int64_t)nside;
  int64_t expected_size = nside_long * nside_long * 12;

  MapData M;
  M.Allocate();
  M.Set_MPIType();
  M.npixels_tot = expected_size;
  vector<map_properties> maps_send; // initialize so we only resize once and keep the memory 
  vector<map_properties> maps_recv;
  vector<float> out_float;
  out_float.reserve(expected_size/commRanks+1); // reserve the expected amount of output size


  string start_file;
  //CHECK 1:  start with a check that all the steps ran and gave an output
  for (int jj=step_start;jj<step_end+1;jj++){
    MPI_Barrier(MPI_COMM_WORLD);
    char step[10*sizeof(char)];
    char mpiioName[512];
    sprintf(step,"%d",jj);
    if (folders){
      strcpy(mpiioName, mpiioName_base);
      strcat(mpiioName, "/step_");
      strcat(mpiioName, step);
      strcat(mpiioName, "/");
      strcat(mpiioName, mpiioName_file);
      strcat(mpiioName, step);
      strcat(mpiioName, ".gio");
    }
    else{
      strcpy(mpiioName, mpiioName_base);
      strcat(mpiioName, "/");
      strcat(mpiioName, mpiioName_file);
      strcat(mpiioName, step);
      strcat(mpiioName, ".gio");
    }
    /*if (check_file(mpiioName)==0){
      if (commRank==0){
        fprintf(stderr,"File does not exist for step %s \n", step);
        fprintf(stderr,"output of check file %d \n", check_file(mpiioName));
        fprintf(stderr,"Filename is %s \n", mpiioName);
      }
      exit(-1);
    }*/
    //assert(check_file(mpiioName)); // assert that the check worked
    //if (jj==step_start){
    //  start_file = mpiioName;
    //}
   
    // if file exists create an hdf5 file for the step  
    //stringstream ss;
    //ss<< commRank;
    //string str_pix = ss.str();
    /*if (commRank==0){
      cout << outfile + "/GO_maps_step_" + step +".hdf5" << "\n";
      string file_name_hdf5 = outfile + "/GO_maps_step_" + step +".hdf5";
      H5File* file = new H5File(file_name_hdf5.c_str(), H5F_ACC_TRUNC);
      delete file;
    }*/

    int status;
    // copy data into map (first step will resize, others should be initialized already)
    status = read_and_redistribute(mpiioName, commRanks, &M, maps_send, maps_recv);
    status = write_files_slab(M, commRank, step, out_float);

    } // close step loop

    M.Deallocate();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

} // close main

