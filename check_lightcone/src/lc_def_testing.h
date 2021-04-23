//#ifndef COSMO_HALO_DEF_H
//#define COSMO_HALO_DEF_H

#include <string>

// List all the halo float fields here
#define N_LC_FLOATS 8

// these are referenced in tree building -- all others are carried allong
enum named_fields {fof_halo_mass=0, fof_halo_center_x=2, fof_halo_center_y=3, fof_halo_center_z=4};

const std::string float_var_names_test[N_LC_FLOATS] = {
	"x",
	"y",
	"z",
	"vx",
	"vy",
	"vz",
	"a",
	"phi"
};

