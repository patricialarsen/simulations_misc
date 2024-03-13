#include <string>
#include <vector>
#include "BasicDefinition.h"


#ifndef HYBRID
#include "halo_def.h"
#else
#include "halo_def.h"
#endif

using namespace std;


class PLHalos {

public:

  int step_number;
  size_t nparticles;
  bool is_allocated;

  vector<int>* pix_index;
  vector<vector<float>* > float_data;
  vector<vector<double>* > double_data;
  vector<vector<int>* > int_data;
  vector<vector<int64_t>* > int64_data;
  vector<vector<uint16_t>* > mask_data;

  vector<int>* rank;
  PLHalos():  nparticles(0), is_allocated(false),\
       step_number(-1)
  {

    float_data.resize(N_FLOATS_HALOS);
    int_data.resize(N_INTS_HALOS);
    int64_data.resize(N_INT64S_HALOS);
    double_data.resize(N_DOUBLES_HALOS);
    mask_data.resize(N_MASKS_HALOS);
  };

  ~PLHalos() { };

  void Allocate(size_t n=0);
  void Deallocate();
  void Resize(size_t);
};

