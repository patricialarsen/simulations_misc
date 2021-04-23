//#ifndef COSMO_HALO_DEF_H
//#define COSMO_HALO_DEF_H

#include <string>

// List all the halo float fields here
#define N_HALO_FLOATS 51

// these are referenced in tree building -- all others are carried allong
enum named_fields {fof_halo_mass=0, fof_halo_center_x=2, fof_halo_center_y=3, fof_halo_center_z=4};

const std::string float_var_names_test[N_HALO_FLOATS] = {
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
	"fof_halo_mean_z",
	"fof_halo_com_x",
	"fof_halo_com_y",
	"fof_halo_com_z",
	"fof_halo_mean_vx",
	"fof_halo_mean_vy",
	"fof_halo_mean_vz",
	"fof_halo_mean_vx",
        "fof_halo_mean_vy",
        "fof_halo_mean_vz",
	"fof_halo_vel_disp",
	"fof_halo_1D_vel_disp",
        "sod_halo_radius",
	"sod_halo_mass",
	"sod_halo_ke",
	"sod_halo_1d_vdisp",
	"sod_halo_max_cir_vel",
	"sod_halo_min_pot_x",
	"sod_halo_min_pot_y",
	"sod_halo_min_pot_z",
	"sod_halo_mean_x",
	"sod_halo_mean_y",
	"sod_halo_mean_z",
	"sod_halo_angmom_x",
	"sod_halo_angmom_y",
	"sod_halo_angmom_z",
	"sod_halo_com_x",
	"sod_halo_com_y",
	"sod_halo_com_z",
	"sod_halo_mean_vx",
	"sod_halo_mean_vy",
	"sod_halo_mean_vz",
	"sod_halo_mean_vx",
        "sod_halo_mean_vy",
        "sod_halo_mean_vz",
	"sod_halo_vel_disp",
	"sod_halo_cdelta",
	"sod_halo_cdelta_error",
	"sod_halo_c_acc_mass",
	"sod_halo_c_peak_mass"
};


const std::string float_var_names_test2[N_HALO_FLOATS] = {
        "fof_halo_mass",
        "fof_halo_ke",
        "fof_halo_center_x",
	"fof_halo_center_y",
        "fof_halo_center_z",
        "fof_halo_angmom_x",
        "fof_halo_angmom_y",
        "fof_halo_angmom_z",//also enter option
        "fof_halo_max_cir_vel",
        "fof_halo_mean_x",
        "fof_halo_mean_y",
        "fof_halo_mean_z",
        "fof_halo_com_x",
        "fof_halo_com_y",
        "fof_halo_com_z",
        "fof_halo_mean_vx", // also com_vx
        "fof_halo_mean_vy",
        "fof_halo_mean_vz",
	"fof_halo_com_vx",
	"fof_halo_com_vy",
        "fof_halo_com_vz",
        "fof_halo_1D_vel_disp", // note repeat, can do this to compare variables
        "fof_halo_1D_vel_disp",
        "sod_halo_radius",
        "sod_halo_mass",
        "sod_halo_ke",
        "sod_halo_1D_vel_disp",
        "sod_halo_max_cir_vel",
        "sod_halo_center_x",
        "sod_halo_center_y",
        "sod_halo_center_z",
        "sod_halo_mean_x",
        "sod_halo_mean_y",
        "sod_halo_mean_z",
	"sod_halo_angmom_x",
	"sod_halo_angmom_y",
	"sod_halo_angmom_z",
        "sod_halo_com_x",
        "sod_halo_com_y",
        "sod_halo_com_z",
        "sod_halo_mean_vx",//also have com_vx
        "sod_halo_mean_vy",
        "sod_halo_mean_vz",
	"sod_halo_com_vx",
	"sod_halo_com_vy",
	"sod_halo_com_vz",
        "sod_halo_1D_vel_disp",
        "sod_halo_cdelta",
        "sod_halo_cdelta_error",
        "sod_halo_c_acc_mass",
        "sod_halo_c_peak_mass"
};

//#endif
