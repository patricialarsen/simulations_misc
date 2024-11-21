//#ifndef HALOS_H
//#define HALOS_H

/*
 * Use this file to access and allocate halo data.
 * Used by merger tree and core tracking code.
 */

#include <string>
#include <vector>
#include <map>
//#include "AllocatedVector.h"
#include "map_def.h"
#include "GenericIO.h"

using namespace std;

// used to carry around the properties for an individual halo


struct map_properties {
  int rank;
  double map_data[N_MAPS];
  int64_t pix_index;
};

struct map_properties_hydro {
  int rank;
  double map_data[N_MAPS_HYDRO];
  int64_t pix_index;
};

class MapDataHydro {

public:

  int step_number;
  size_t npixels;
  size_t npixels_tot;
  MPI_Datatype map_properties_hydro_MPI_Type;
  bool is_allocated;

  vector<int64_t> pix_index;
  vector<vector<double>* > double_data;

  // destination rank for redistribution

  MapDataHydro():  npixels(0), is_allocated(false),\
       step_number(-1), npixels_tot(0)
  {
    double_data.resize(N_MAPS_HYDRO);
  };

  ~MapDataHydro() { };

  // check if others are needed 
  void Allocate(size_t n=0);
  void Deallocate();
  void Resize(size_t);
  void Clear(size_t);
  void Set_MPIType();
  map_properties_hydro GetProperties(size_t);
  void Assign(map_properties_hydro, size_t);

};

class MapData {

public:

  int step_number;
  size_t npixels;
  size_t npixels_tot;

  MPI_Datatype map_properties_MPI_Type;
  bool is_allocated;

  vector<int64_t> pix_index;

  vector<vector<double>* > double_data;


  // destination rank for redistribution
  //vector<int> rank;

  MapData():  npixels(0), is_allocated(false),\
       step_number(-1), npixels_tot(0)
  {
    double_data.resize(N_MAPS);
  };

  ~MapData() { };

  // check if others are needed
  void Allocate(size_t n=0);
  void Deallocate();
  void Resize(size_t);
  void Clear(size_t);
  void Set_MPIType();
  map_properties GetProperties(size_t);
  void Assign(map_properties, size_t);

};


