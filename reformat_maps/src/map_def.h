#include <string>
#include <vector>

#define N_MAPS 10
//#define N_MAPS 3

#define N_MAPS_HYDRO 10
using namespace std;

const string map_names_hydro[3] = {
        "phi",
	"vel",
	"rho",
};

//const string map_names_hydro[N_MAPS_HYDRO] = {
const string map_names_go[N_MAPS] = {
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

