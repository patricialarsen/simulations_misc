#!/bin/bash/python

import sys
sys.path.append('/home/prlarsen/usr/genericio/python')
import genericio as gio
import numpy as np
import matplotlib.pyplot as plt

file_old = '/projects/DarkUniverse_esp/prlarsen/hacc_tests/center_finder/analysis_final_tmp/m000p-499.haloproperties'
file_new = '/projects/DarkUniverse_esp/prlarsen/hacc_tests/center_finder/analysis_final3/m000p-499.haloproperties'

def return_idx():
    tag_new = gio.gio_read(file_new,'fof_halo_tag')
    tag_old = gio.gio_read(file_old,'fof_halo_tag')
    idx_new = np.argsort(tag_new)
    idx_old = np.argsort(tag_old)
    print((tag_new[idx_new]==tag_old[idx_old]).all())
    return idx_new, idx_old

idx_new,idx_old = return_idx()

def check_prop(var,idx_new,idx_old):
    var_new = gio.gio_read(file_new,var)
    var_old = gio.gio_read(file_old,var)
    ident = (var_new[idx_new]==var_old[idx_old]).all()
    if ident:
        print("values are identical for variable ",var)
    else:
        print("values are not identical for variable ",var)
        print(np.sum(var_new[idx_new]!=var_old[idx_old]))
        print(np.max(var_new[idx_new]-var_old[idx_old]))
        print(np.min(var_new[idx_new]-var_old[idx_old]))
        print(np.max((var_new[idx_new]-var_old[idx_old])/var_old[idx_old]))
        print(np.min((var_new[idx_new]-var_old[idx_old])/var_old[idx_old]))
    return 

  
prop_list1 = ['fof_halo_count','fof_halo_tag','fof_halo_mass','fof_halo_ke','fof_halo_center_x','fof_halo_center_y','fof_halo_center_z','fof_halo_angmom_x','fof_halo_angmom_y','fof_halo_angmom_z']


prop_list2 = ['fof_halo_max_cir_vel','fof_halo_com_x','fof_halo_com_y','fof_halo_com_z','fof_halo_mean_vx','fof_halo_mean_vy','fof_halo_mean_vz','fof_halo_vel_disp','fof_halo_1D_vel_disp']


prop_list3=[]
for i in ['S','R']:
    for j in ['1','2','3']:
        for k in ['X','Y','Z']:
            prop_list3.append('fof_halo_eig'+i+j+k)

prop_list4= ['sod_halo_count','sod_halo_radius','sod_halo_mass','sod_halo_ke','sod_halo_1d_vdisp','sod_halo_max_cir_vel','sod_halo_min_pot_x','sod_halo_min_pot_y','sod_halo_min_pot_z','sod_halo_angmom_x','sod_halo_angmom_y','sod_halo_angmom_z']

prop_list5 = ['sod_halo_mean_x','sod_halo_mean_y','sod_halo_mean_z','sod_halo_mean_vx','sod_halo_mean_vy','sod_halo_mean_vz','sod_halo_vel_disp','sod_halo_cdelta','sod_halo_cdelta_error','sod_halo_c_acc_mass','sod_halo_c_peak_mass']

prop_list6=[]
for i in ['S','R']:
    for j in ['1','2','3']:
        for k in ['X','Y','Z']:
            prop_list6.append('sod_halo_eig'+i+j+k)

