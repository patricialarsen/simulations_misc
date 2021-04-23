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
#include "lc_def_testing.h"
#include "GenericIO.h"

using namespace std;

// used to carry around the properties for an individual lc particle
struct lc_properties_test {
  int64_t id;
  int32_t replication;
  int32_t rank;
  float float_data[N_LC_FLOATS];
};

class LC_test {

public:

  bool is_allocated;
  int step_number;
  size_t num_parts;
  float particle_mass;
  MPI_Datatype LC_properties_MPI_Type;

  vector<int64_t>* id;
  vector<int32_t>* replication;
  vector<vector<float>* > float_data;

  // destination rank for redistribution
  vector<int32_t>* rank; 

  // map for the halo tag to the index in the local arrays
  map<int64_t,int>* tag2idx;

  // halo particle membership
  //char*  buffer;                            // a memory buffer for the pids and tags
 // size_t buffer_size;                       // size (in bytes) of the buffer
 // allocated_vector<int64_t>* tag;           // ID for halos on this rank,
 // allocated_vector<int64_t>* pid;           // ID for particles on this rank

  LC_test() :  num_parts(0), particle_mass(0.f), is_allocated(false),\
      step_number(-1)
  {
    float_data.resize(N_LC_FLOATS);
  };

  ~LC_test() { };

  void Allocate(size_t n=0);
  void Deallocate();
  lc_properties_test GetProperties(size_t idx);
  void PushBack(lc_properties_test);
  void Resize(size_t);
  void Erase(size_t);
  void Set_MPIType();
};

//#endif
