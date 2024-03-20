#include <string>
#define N_FLOATS 15
#define N_MASKS 1
#define N_INTS 1 
#define N_INT64S 2
#define N_DOUBLES 0 

enum named_fields {x=0, y=1, z=2, a=7, mask=0};

// first three have to be x,y,z 
// next three are vx,vy,vz
const std::string float_names[N_FLOATS] = {
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
	"phi",
        "a",
        "dt",
        "mass",
        "uu",
        "mu",
        "rho",
        "zmet",
        "yhe"
};

const std::string int_names[N_INTS] = {
       "replication"
};

const std::string int64_names[N_INT64S] = {
       "id",
       "fof_halo_tag"
};

const std::string mask_names[N_MASKS] = {
       "mask"
};

const std::string double_names[N_DOUBLES] = {
};
