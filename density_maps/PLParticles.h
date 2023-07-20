#include <string>
#include <vector>
#include "BasicDefinition.h"


#ifndef HYBRID
#include "particle_def.h"
#else
#include "particle_def_hydro.h"
#endif

using namespace std;


class PLParticles {

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
  PLParticles():  nparticles(0), is_allocated(false),\
       step_number(-1)
  {

    float_data.resize(N_FLOATS);
    int_data.resize(N_INTS);
    int64_data.resize(N_INT64S);
    double_data.resize(N_DOUBLES);
    mask_data.resize(N_MASKS);
  };

  ~PLParticles() { };

  void Allocate(size_t n=0);
  void Deallocate();
  void Resize(size_t);
};

