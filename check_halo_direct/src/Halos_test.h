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
};

class Halos_test {

public:

  bool has_sod;
  bool is_allocated;
  int step_number;
  int descendant_step_number;
  size_t num_halos;
  float particle_mass;

  vector<int64_t>* fof_halo_tag;
  vector<int32_t>* fof_halo_count;
  vector<int64_t>* sod_halo_count;

  vector<vector<float>* > float_data;

  // destination rank for redistribution
  vector<int32_t>* rank; 


  // map for the halo tag to the index in the local arrays
  map<int64_t,int>* tag2idx;

  // halo particle membership
  char*  buffer;                            // a memory buffer for the pids and tags
  size_t buffer_size;                       // size (in bytes) of the buffer
  allocated_vector<int64_t>* tag;           // ID for halos on this rank,
  allocated_vector<int64_t>* pid;           // ID for particles on this rank

  Halos_test() : has_sod(true), num_halos(0), particle_mass(0.f), is_allocated(false),\
      buffer(NULL), buffer_size(0), step_number(-1),\
      descendant_step_number(-1) 
  {
    float_data.resize(N_HALO_FLOATS);
  };

  ~Halos_test() { };

  void Allocate(size_t n=0);
  void Deallocate();
  halo_properties_test GetProperties(size_t idx);
  void PushBack(halo_properties_test);
  void Resize(size_t);
  void Erase(size_t);
};

//#endif
