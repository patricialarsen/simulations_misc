/*
 * Use this file to access and allocate SZ-related data.
 */

#include <string>
#include <vector>
#include <map>
#include "AllocatedVector.h"
#include "GenericIO.h"

using namespace std;

struct sz_props{
    int64_t pix_num;
    double tsz;
    double ksz;
    int rank;
};

class SZ_class {

public:

  bool is_allocated;
  int step_number;
  size_t num_parts;
  MPI_Datatype sz_properties_MPI_Type;

  vector<int64_t>* pix_num;
  vector<double>* tsz; 
  vector<double>* ksz;
  // destination rank for redistribution
  vector<int>* rank; 

  SZ_class() :  num_parts(0), is_allocated(false),\
      step_number(-1);
  ~SZ_class() { };

  void Allocate(size_t n=0);
  void Deallocate();
  sz_props GetProperties(size_t idx);
  void Resize(size_t);
  void Set_MPIType();
};

