//#ifndef HALOS_H
//#define HALOS_H

/*
 * Use this file to access and allocate halo data.
 * Used by merger tree and core tracking code.
 */

#include <string>
#include <vector>
#include <map>
#include "AllocatedVector.h"
#include "halo_def_LJ.h"
#include "GenericIO.h"

using namespace std;

// used to carry around the properties for an individual halo

#define N_HALO_FLOATS 32
#define N_LIGHTCONE_FLOATS 7
#define N_HALO_FLOATS_E 18


struct halo_properties {
  int64_t fof_halo_tag;
  int rank;
  int pix_index;
  float lightcone_data[N_LIGHTCONE_FLOATS];
  float float_data[N_HALO_FLOATS];
  float ellipticity_data[N_HALO_FLOATS_E];
};


class Halos_test {

public:

  int step_number;
  size_t num_halos;
  float particle_mass;
  MPI_Datatype halo_properties_MPI_Type;
  bool is_allocated;


  vector<int>* pix_index;
  vector<int64_t>* fof_halo_tag;
  vector<vector<float>* > lightcone_data;
  vector<vector<float>* > float_data;
  vector<vector<float>* > ellipticity_data;

  // destination rank for redistribution
  vector<int>* rank; 

  // map for the halo tag to the index in the local arrays
  //map<int64_t,int>* tag2idx;

  Halos_test():  num_halos(0), particle_mass(0.f), is_allocated(false),\
       step_number(-1)
  {
    float_data.resize(N_HALO_FLOATS);
    ellipticity_data.resize(N_HALO_FLOATS_E);
    lightcone_data.resize(N_LIGHTCONE_FLOATS);
  };

  ~Halos_test() { };

  void Allocate(size_t n=0);
  void Deallocate();
  halo_properties GetProperties(size_t idx);
  void PushBack(halo_properties);
  void Resize(size_t);
  void Assign(halo_properties, size_t);
  void Erase(size_t);
  void Set_MPIType();
};


