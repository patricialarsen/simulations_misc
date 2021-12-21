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
#include "halo_def_testing.h"
#include "GenericIO.h"

using namespace std;

// used to carry around the properties for an individual halo
struct halo_properties_test {
  int64_t fof_halo_tag;
  int32_t fof_halo_count;
  int64_t sod_halo_count;
  int32_t rank;
  float float_data[N_HALO_FLOATS];
  float ellipticity_data[N_HALO_FLOATS_E];
};

struct sod_binproperties_test {
  int64_t fof_halo_bin_tag;
  int32_t sod_halo_bin;
  float float_data[N_HALO_FLOATS_SOD];
  int int_data[N_HALO_INTS_SOD];
  int32_t rank;
};

struct particles_test {
  int64_t fof_halo_tag;
  int64_t id;
  int32_t rank;
};

class Halos_test {

public:

  bool has_sod;
  bool is_allocated;
  int step_number;
  size_t num_halos;
  float particle_mass;
  MPI_Datatype halo_properties_MPI_Type;

  vector<int64_t>* fof_halo_tag;
  vector<int32_t>* fof_halo_count;
  vector<int64_t>* sod_halo_count;

  vector<vector<float>* > float_data;
  vector<vector<float>* > ellipticity_data;

  // destination rank for redistribution
  vector<int32_t>* rank; 

  // map for the halo tag to the index in the local arrays
  map<int64_t,int>* tag2idx;

  Halos_test() : has_sod(true), num_halos(0), particle_mass(0.f), is_allocated(false),\
       step_number(-1)
  {
    float_data.resize(N_HALO_FLOATS);
    ellipticity_data.resize(N_HALO_FLOATS_E);
  };

  ~Halos_test() { };

  void Allocate(size_t n=0);
  void Deallocate();
  halo_properties_test GetProperties(size_t idx);
  void PushBack(halo_properties_test);
  void Resize(size_t);
  void Erase(size_t);
  void Set_MPIType();
};


class SODBins_test {

public:

  bool is_allocated;
  int step_number;
  size_t num_halos;
  float particle_mass;
  MPI_Datatype sod_binproperties_MPI_Type;

  vector<int64_t>* fof_halo_bin_tag;
  vector<int32_t>* sod_halo_bin;

  vector<vector<float>* > float_data;
  vector<vector<int>* > int_data;
  vector<int32_t>* rank;

  map<int64_t,int>* tag2idx;

  SODBins_test() : num_halos(0), particle_mass(0.f), is_allocated(false),\
      step_number(-1)
  {
    float_data.resize(N_HALO_FLOATS_SOD);
    int_data.resize(N_HALO_INTS_SOD);
  };

  ~SODBins_test() { };

  void Allocate(size_t n=0);
  void Deallocate();
  sod_binproperties_test GetProperties(size_t idx);
  void PushBack(sod_binproperties_test);
  void Resize(size_t);
  void Erase(size_t);
  void Set_MPIType();
};


class Particles_test {

public:

  bool is_allocated;
  int step_number;
  size_t num_halos;
  MPI_Datatype particles_test_MPI_Type;

  vector<int64_t>* fof_halo_tag;
  vector<int64_t>* id;

  vector<int32_t>* rank;

  map<int64_t,int>* tag2idx;

  Particles_test() : num_halos(0), is_allocated(false),\
      step_number(-1) {};

  ~Particles_test() { };

  void Allocate(size_t n=0);
  void Deallocate();
  particles_test GetProperties(size_t idx);
  void PushBack(particles_test);
  void Resize(size_t);
  void Erase(size_t);
  void Set_MPIType();
};



//#endif
