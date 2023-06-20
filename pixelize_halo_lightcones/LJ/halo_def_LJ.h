#include <string>

// List all the halo float fields here
#define N_HALO_FLOATS 32
#define N_LIGHTCONE_FLOATS 7
#define N_HALO_FLOATS_E 18

// these are referenced in tree building -- all others are carried along
// must retain this order for these codes - do not change this line
enum named_fields {fof_halo_mass=0, fof_halo_center_x=2, fof_halo_center_y=3, fof_halo_center_z=4};


const std::string lightcone_var_names[N_LIGHTCONE_FLOATS] = {
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "a"
};


const std::string float_var_names[N_HALO_FLOATS] = {
	"fof_halo_mass",
	"fof_halo_ke",
	"fof_halo_center_x",
	"fof_halo_center_y",
	"fof_halo_center_z",
	"fof_halo_angmom_x",
	"fof_halo_angmom_y",
	"fof_halo_angmom_z", 
	"fof_halo_max_cir_vel",
	"fof_halo_mean_x",
	"fof_halo_mean_y",
	"fof_halo_mean_z" ,
	"fof_halo_mean_vx",
	"fof_halo_mean_vy",
	"fof_halo_mean_vz",
	"fof_halo_vel_disp",
        "sod_halo_radius",
	"sod_halo_mass",
	"sod_halo_ke",
	"sod_halo_vel_disp",
	"sod_halo_max_cir_vel",
	"sod_halo_mean_x",
	"sod_halo_mean_y",
	"sod_halo_mean_z",
	"sod_halo_angmom_x",
	"sod_halo_angmom_y",
	"sod_halo_angmom_z",
	"sod_halo_mean_vx",
	"sod_halo_mean_vy",
	"sod_halo_mean_vz",
	"sod_halo_cdelta",
	"sod_halo_cdelta_error",
};



const std::string float_var_names_ellipticity[N_HALO_FLOATS_E] = {
        "fof_halo_eigS1X",
        "fof_halo_eigS1Y",
        "fof_halo_eigS1Z",
        "fof_halo_eigS2X",
        "fof_halo_eigS2Y",
        "fof_halo_eigS2Z",
        "fof_halo_eigS3X",
        "fof_halo_eigS3Y",
        "fof_halo_eigS3Z",
        "sod_halo_eigS1X",
        "sod_halo_eigS1Y",
        "sod_halo_eigS1Z",
        "sod_halo_eigS2X",
        "sod_halo_eigS2Y",
        "sod_halo_eigS2Z",
        "sod_halo_eigS3X",
        "sod_halo_eigS3Y",
        "sod_halo_eigS3Z"
};


