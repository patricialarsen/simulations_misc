#include <string>
#include <vector>

#define N_MAPS 3
#define N_MAPS_HYDRO 10
using namespace std;

const string map_names_go[N_MAPS] = {
        "phi",
	"vel",
	"rho",
};

const string map_names_hydro[N_MAPS_HYDRO] = {
	"phi",
	"vel",
	"rho",
	"xray1",
	"xray2",
	"xray3",
	"xray4",
	"temp",
	"ksz",
	"tsz",
};


class MapData {

public:

  int step_number;
  size_t npixels;
  bool is_allocated;

  // vector of pixel indices
  vector<int64_t> pix_index;
  // vector of vector pointers
  vector<vector<double>* > double_data;

  MapData():  npixels(0), is_allocated(false),\
       step_number(-1)
  {
    // resizing the vector of pointers to the number of maps
    double_data.resize(N_MAPS);
  };

  ~MapData() { };

  void Allocate(size_t n=0);
  void Deallocate();
  void Resize(size_t);
  void Clear(size_t);

};

class MapDataHydro {

public:

  int step_number;
  size_t npixels;
  bool is_allocated;

  vector<int64_t> pix_index;
  vector<vector<double>* > double_data;

  MapDataHydro():  npixels(0), is_allocated(false),\
       step_number(-1)
  {
    double_data.resize(N_MAPS_HYDRO);
  };

  ~MapDataHydro() { };

  void Allocate(size_t n=0);
  void Deallocate();
  void Resize(size_t);
  void Clear(size_t);

};

