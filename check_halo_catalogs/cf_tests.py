#!/bin/bash/python

import sys
sys.path.append('/home/prlarsen/usr/genericio/python')
import genericio as gio
import numpy as np
import matplotlib.pyplot as plt

file_old = '/gpfs/mira-fs0/projects/DarkUniverse_esp/prlarsen/hacc_tests/center_finder/analysis_trunk/m000p-499.haloproperties'
file_new = '/gpfs/mira-fs0/projects/DarkUniverse_esp/prlarsen/hacc_tests/center_finder/analysis/m000p-499.haloproperties'

def check_tags():
    ''' Double check ordering is the same from the tags '''
    tag_old = gio.gio_read(file_old,'fof_halo_tag')
    tag_new = gio.gio_read(file_new,'fof_halo_tag')
    if (tag_old==tag_new).all():
         print("ordering is the same")
    else:
         print("ordering has changed - double check this before continuing")
    return 

def make_mask():
    ''' Make a 10^12 mass cut mask from the FOFs '''
    fof_old = gio.gio_read(file_old,'fof_halo_mass')
    return (fof_old>1.e12), fof_old[fof_old>1.e12]

def read_mass_cdelta():
    ''' read and mask the SOD mass and concentration '''
    mask12, fof_old = make_mask()
    print("made mask")
    sod_old = gio.gio_read(file_old,'sod_halo_mass')[mask12]
    sod_new = gio.gio_read(file_new,'sod_halo_mass')[mask12]
    cdelta_old = gio.gio_read(file_old,'sod_halo_cdelta')[mask12]
    cdelta_new = gio.gio_read(file_new,'sod_halo_cdelta')[mask12]
    return fof_old, sod_old, sod_new, cdelta_old, cdelta_new

def read_centers():
    mask12, fof_old = make_mask()
    x_old = gio.gio_read(file_old,'fof_halo_center_x')[mask12]
    y_old = gio.gio_read(file_old,'fof_halo_center_y')[mask12]
    z_old = gio.gio_read(file_old,'fof_halo_center_z')[mask12]
    x_new = gio.gio_read(file_new,'fof_halo_center_x')[mask12]
    y_new = gio.gio_read(file_new,'fof_halo_center_y')[mask12]
    z_new = gio.gio_read(file_new,'fof_halo_center_z')[mask12]
    radius = gio.gio_read(file_old,'sod_halo_radius')[mask12]
    return fof_old, x_old,y_old,z_old,x_new,y_new,z_new, radius

def make_dists():
    fof_old, x_old,y_old,z_old,x_new,y_new,z_new, radius = read_centers()
    dist = (x_new-x_old)**2 + (y_new -y_old)**2 + (z_new-z_old)**2
    dist = np.sqrt(dist)
    return fof_old, dist , radius

def make_dist_plots(fof_old,dist,radius,mass_cut):
    mask_final = (fof_old>mass_cut)
    x_range  = np.linspace(np.min(np.log10(fof_old[mask_final])),np.max(np.log10(fof_old[mask_final])),100)
    plt.figure()
    plt.scatter(np.log10(fof_old[mask_final]),(dist/radius)[mask_final],s=0.1)
    plt.xlabel('log10(M_FOF)')
    plt.ylabel('shift distance/SO radius')
    plt.savefig('fig_'+str(mass_cut)+'_distance_full.png',dpi=200)
    plt.show()
    plt.figure()
    plt.scatter(np.log10(fof_old[mask_final]),(dist/radius)[mask_final],s=0.1)
    plt.xlabel('log10(M_FOF)')
    plt.ylabel('shift distance/SO radius')
    plt.ylim([-0.1,0.1])
    plt.fill_between(x_range, -0.01,0.01,alpha=0.1,color='k')
    plt.savefig('fig_'+str(mass_cut)+'_distance_zoom.png',dpi=200)
    plt.show()
    return 


def make_sod_plots(fof_old,sod_old,sod_new,cdelta_old,cdelta_new, mass_cut):
    ''' Make plots related to sod mass and cdelta '''
    mask_final = (fof_old>mass_cut)
    x_range  = np.linspace(np.min(np.log10(fof_old[mask_final])),np.max(np.log10(fof_old[mask_final])),100)
    ratio_mass = (sod_new-sod_old)/sod_old
    ratio_cdelta = (cdelta_new - cdelta_old)/cdelta_old
    plt.figure()
    plt.scatter(np.log10(sod_old[mask_final]),np.log10(sod_new[mask_final]),s=0.1)
    plt.ylabel('M200 (global potential in Msun/h)')
    plt.xlabel('M200 (most bound particle in Msun/h)')
    plt.savefig('fig_'+str(mass_cut)+'_M200.png',dpi=200)
    plt.show()
    plt.figure()
    plt.scatter(cdelta_old[mask_final], cdelta_new[mask_final],s=0.1)
    plt.xlabel('concentration (most bound particle)')
    plt.ylabel('concentration (global potential)')
    plt.savefig('fig_'+str(mass_cut)+'_cdelta.png',dpi=200)
    plt.show()
    plt.figure()
    plt.scatter(np.log10(fof_old[mask_final]),ratio_mass[mask_final],s=0.1)
    plt.ylabel('(M200 global - M200 MBP)/(M200 MBP)')
    plt.xlabel('log10(M_FOF)')
    plt.savefig('fig_'+str(mass_cut)+'_M200_ratio_full.png',dpi=200)
    plt.show()
    plt.figure()
    plt.scatter(np.log10(fof_old[mask_final]),ratio_mass[mask_final],s=0.1)
    plt.ylim([-0.1,0.1])
    plt.fill_between(x_range, -0.01,0.01,alpha=0.1,color='k')
    plt.ylabel('(M200 global - M200 MBP)/(M200 MBP)')
    plt.xlabel('log10(M_FOF)')
    plt.savefig('fig_'+str(mass_cut)+'_M200_ratio_zoom.png',dpi=200)
    plt.show()
    plt.figure()
    plt.scatter(np.log10(fof_old[mask_final]),ratio_cdelta[mask_final],s=0.1)
    plt.ylabel('(cdelta global - cdelta MBP)/(cdelta MBP)')
    plt.xlabel('log10(M_FOF)')
    plt.savefig('fig_'+str(mass_cut)+'_cdelta_ratio_zoom.png',dpi=200)
    plt.show()
    plt.figure()
    plt.scatter(np.log10(fof_old[mask_final]),ratio_cdelta[mask_final],s=0.1)
    plt.fill_between(x_range, -0.01,0.01,alpha=0.1,color='k')
    plt.ylim([-0.1,0.1])
    plt.ylabel('(cdelta global - cdelta MBP)/(cdelta MBP)')
    plt.xlabel('log10(M_FOF)')
    plt.savefig('fig_'+str(mass_cut)+'_cdelta_ratio_zoom.png',dpi=200)
    plt.show()
    return 

def percentages(fof_old,sod_old,sod_new,cdelta_old,cdelta_new,dist,radius, mass_cut):
    mask_final = (fof_old>mass_cut)
    ratio_mass = (sod_new-sod_old)/sod_old
    ratio_cdelta = (cdelta_new - cdelta_old)/cdelta_old
    ratio_dist = (dist/radius)
    perc_mass = np.sum(ratio_mass[mask_final]<0.01)/float(len(ratio_mass[mask_final]))*100.
    perc_cdelta = np.sum(ratio_cdelta[mask_final]<0.01)/float(len(ratio_cdelta[mask_final]))*100.
    perc_dist = np.sum(ratio_dist[mask_final]<0.01)/float(len(ratio_dist[mask_final]))*100.
    print("percentage of objects with ratio_mass less than 1% = ", perc_mass)
    print("percentage of objects with ratio_cdelta less than 1% = ", perc_cdelta)
    print("percentage of objects with ratio_dist less than 1% = ", perc_dist)
    perc_mass = np.sum(ratio_mass[mask_final]==0.0)/float(len(ratio_mass[mask_final]))*100.
    perc_cdelta = np.sum(ratio_cdelta[mask_final]==0.0)/float(len(ratio_cdelta[mask_final]))*100.
    perc_dist = np.sum(ratio_dist[mask_final]==0.0)/float(len(ratio_dist[mask_final]))*100.
    print("percentage of objects with ratio_mass of 0 = ", perc_mass)
    print("percentage of objects with ratio_cdelta of 0 = ", perc_cdelta)
    print("percentage of objects with ratio_dist of 0 = ", perc_dist)

    mean_mass = np.mean(ratio_mass[mask_final])
    mean_cdelta = np.mean(ratio_cdelta[mask_final])
    std_mass = np.std(ratio_mass[mask_final])
    std_cdelta = np.std(ratio_cdelta[mask_final])
    print("mean for mass = ",mean_mass)
    print("std for mass = ",std_mass)
    print("mean for c = ",mean_cdelta)
    print("std for c = ",std_cdelta)

    return perc_mass, perc_cdelta


def above_below(fof_old,sod_old,sod_new,cdelta_old,cdelta_new,dist,radius, mass_cut):
    mask_final = (fof_old>mass_cut)
    ratio_mass = (sod_new-sod_old)/sod_old
    ratio_cdelta = (cdelta_new - cdelta_old)/cdelta_old
    ratio_dist = (dist/radius)
    print("number with sod_mass higher in global potential = ", np.sum(ratio_mass[mask_final]>0))
    print("number with sod_mass lower in global potential = ", np.sum(ratio_mass[mask_final]<0))
    print("number with sod_cdelta higher in global potential = ", np.sum(ratio_cdelta[mask_final]>0))
    print("number with sod_cdelta lower in global potential = ", np.sum(ratio_cdelta[mask_final]<0))
    return

