#include <string>
#define N_FLOATS_HALOS 5
#define N_MASKS_HALOS 0
#define N_INTS_HALOS 0 
#define N_INT64S_HALOS 1 
#define N_DOUBLES_HALOS 0 

//enum named_fields_halos {x=0, y=1, z=2, a=7, mask=0};

// first three have to be x,y,z 
// next three are vx,vy,vz
const std::string float_names_halos[N_FLOATS_HALOS] = {
	"fof_halo_center_x",
	"fof_halo_center_y",
        "fof_halo_center_z",
        "sod_halo_mass",
	"sod_halo_radius" 
};
//sod_halo_R500c 
//sod_halo_radius - r200c
//sod_halo_R200m 

const std::string int_names_halos[N_INTS_HALOS] = {
};

const std::string int64_names_halos[N_INT64S_HALOS] = {
       "fof_halo_tag"
};

const std::string mask_names_halos[N_MASKS_HALOS] = {
};

const std::string double_names_halos[N_DOUBLES_HALOS] = {
};

