#include <cstdlib>
#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdexcept>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <algorithm>    // std::sort
#include <sstream>
#include <omp.h>
#include <sys/stat.h>

// Generic IO
#include "GenericIO.h"
#include "Partition.h"

#include "Halos_test.h"

//#include "MurmurHashNeutral2.cpp" 
#include "routines.h"

/// Assumes ids and  are consistent but ordering may have changed
/// Redistributes aongst ranks and sorts the values for a one-to-one check of the changes in 
/// several halo outputs


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

int main( int argc, char** argv ) {

  MPI_Init( &argc, &argv );
  Partition::initialize();
  GenericIO::setNaturalDefaultPartition();

  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  string path1, path2;
  int halo_opt, sod_opt, part_opt, hydro_opt;
  float lim, mass_min, mass_max, box_size;
  // read params file to give options for what to run 
  if (rank==0){
     fstream params;
     params.open("params.txt", ios::in);
     if (!params) {
        cout << "Parameter file doesn't exist";
     }
     else {
	string::size_type sz;
	string tmp;
	getline(params, path1);
	getline(params, path2);
	getline(params,tmp);
	halo_opt = stoi(tmp,&sz);
	getline(params,tmp);
	getline(params,tmp);
        sod_opt = stoi(tmp,&sz);
	getline(params,tmp);
	getline(params,tmp);
	part_opt = stoi(tmp,&sz);
	getline(params,tmp);
	getline(params,tmp);
	box_size = stof(tmp,&sz);
	getline(params,tmp);
	getline(params,tmp);
	lim = stof(tmp,&sz);
	getline(params,tmp);
	getline(params,tmp);
        mass_min = stof(tmp,&sz);
	getline(params,tmp);
	getline(params,tmp);
	mass_max = stof(tmp,&sz);
        getline(params,tmp);
        getline(params,tmp);
        hydro_opt = stoi(tmp,&sz);
	params.close();

  }
  }
  int path1_size = path1.size();
  int path2_size = path2.size();

  MPI_Bcast(&path1_size, 1, MPI_INT, 0 , Partition::getComm());
  MPI_Bcast(&path2_size, 1, MPI_INT, 0 , Partition::getComm());

  if (rank!=0){
  path1.resize(path1_size);
  path2.resize(path2_size);
  }
  MPI_Bcast(&path1[0], path1.size(), MPI_CHAR, 0,Partition::getComm());
  MPI_Bcast(&path2[0], path2.size(), MPI_CHAR, 0,Partition::getComm());
  MPI_Bcast(&halo_opt, 1, MPI_INT, 0,Partition::getComm());
  MPI_Bcast(&sod_opt, 1, MPI_INT, 0,Partition::getComm());
  MPI_Bcast(&part_opt, 1, MPI_INT, 0,Partition::getComm());
  MPI_Bcast(&box_size, 1, MPI_FLOAT, 0,Partition::getComm());
  MPI_Bcast(&lim, 1, MPI_FLOAT, 0,Partition::getComm());
  MPI_Bcast(&mass_min, 1, MPI_FLOAT, 0,Partition::getComm());
  MPI_Bcast(&mass_max, 1, MPI_FLOAT, 0,Partition::getComm());

  if (rank==0){
  cout << " Parameter list"<< endl;
  cout << " " << path1 << endl;
  cout << " " << path2 << endl;
  cout << " " << halo_opt << endl;
  cout << " " << sod_opt << endl;
  cout << " " << part_opt << endl;
  cout << " " << box_size << endl;
  cout << " " << lim << endl;
  cout << " " << mass_min <<endl;
  cout << " " << mass_max << endl;
  cout << endl;
  }


  string fof_file, fof_file2;
  fof_file = path1+ "/m000p-499.haloproperties";
  fof_file2 = path2 + "/m000p-499.haloproperties";
  string sod_file, sod_file2;
  sod_file = path1+ "/m000p-delta200.0-499.sodpropertybins";
  sod_file2 = path2 + "/m000p-delta200.0-499.sodpropertybins";
  string part_file, part_file2;
  part_file = path1+ "/m000p-499.bighaloparticles#0";
  part_file2 = path2 + "/m000p-499.bighaloparticles#0";

  map<int64_t,int> tag_map_main;
  map<int64_t,int>* tag_map = &tag_map_main; // halos passing the mass thresholds

  int err=0;

  if (halo_opt==0){
    err = perform_halo_check (fof_file, fof_file2, lim, mass_min, mass_max,tag_map);
  }
  else if (halo_opt==1)
    err = match_pos(fof_file, fof_file2, lim, box_size, mass_min, mass_max);
  else if (halo_opt==2)
    err = compare_dist(fof_file, fof_file2, lim);
  if ((sod_opt==0)&&(halo_opt==0))
    err = sodbin_check(sod_file,sod_file2,lim,tag_map);
  if (part_opt==0)
    err = part_check(part_file, part_file2);
  if (part_opt==1){
    err = part_check(part_file, part_file2);
    part_file = path1+ "/m000p-499.haloparticles#0";
    part_file2 = path2 + "/m000p-499.haloparticles#0";
    err = part_check(part_file, part_file2);
    part_file = path1+ "/m000p-499.haloparticletags#0";
    part_file2 = path2 + "/m000p-499.haloparticletags#0";
    err = part_check(part_file, part_file2);
  }


  Partition::finalize();
  MPI_Finalize();
  return 0;
}
