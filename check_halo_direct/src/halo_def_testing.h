//#ifndef COSMO_HALO_DEF_H
//#define COSMO_HALO_DEF_H

#include <string>

// List all the halo float fields here
#define N_HALO_FLOATS 6

// these are referenced in tree building -- all others are carried allong
enum named_fields {fof_halo_mass=0, fof_halo_center_x=2, fof_halo_center_y=3, fof_halo_center_z=4};

const std::string float_var_names_test[N_HALO_FLOATS] = {
  "fof_halo_mass",
  "fof_halo_ke",
  "fof_halo_center_x",
  "fof_halo_center_y",
  "fof_halo_center_z",
  "fof_halo_angmom_x"
};


const std::string float_var_names_test2[N_HALO_FLOATS] = {
  "fof_halo_mass",
  "fof_halo_ke",
  "fof_halo_center_x",
  "fof_halo_center_y",
  "fof_halo_center_z",
  "fof_halo_com_angmom_x"
};

//#endif
