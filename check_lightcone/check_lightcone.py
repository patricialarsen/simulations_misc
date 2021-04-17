#!/bin/python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os 
import sys
sys.path.append('/home/prlarsen/usr/genericio/python')
import genericio as gio



#TODO: change genericio directory and output directory
#os.chdir('/home/prlarsen/usr/genericio/python')
#import genericio as gio
#os.chdir('/home/prlarsen/lightcone_analysis/outputs')


def compute_mask(nside=4096,octant="upper"):
    ''' computing octant mask'''
    x,y,z = hp.pix2vec(nside,np.arange(hp.nside2npix(nside)))
    if (octant == "lower"):
        mask_oct = (x>=0)*(y>=0)*(z<=0)
    else:
        mask_oct = (x>=0)*(y>=0)*(z>=0)
    return mask_oct

def scale_map(map_input,mask_oct,nside=4096):
    mean_input = np.mean(map_input[mask_oct])
    map_scaled = np.zeros(hp.nside2npix(nside))
    map_scaled[mask_oct] = (map_input[mask_oct] - mean_input)/mean_input
    map_scaled[mask_oct==False] = 0
    return map_scaled


def check_file(input_file):
    gio.gio_inspect(input_file)
    return

def read_in_lightcone(input_file,vars=['x','y','z','vx','vy','vz','id']):
    ''' read in lightcone as a dictionary '''
    lightcone_object={}
    for i in vars:
        lightcone_object[i]=gio.gio_read(input_file,i)
    return lightcone_object

def compare_ids(object1,object2):
    ''' search for missing id values '''
    diff_id = set(object1['id']).difference(object2['id'])
    diff_id2 = set(object2['id']).difference(object1['id'])
    return diff_id, diff_id2

def create_mask(object1,diff_id):
    ''' given a set of missing ids, create a mask for these '''
    mask_id = (object1['id']==-1)
    for i in diff_id:
        mask_id += (object1['id']==i)
    return mask_id

def missing_distances(input_file1,input_file2, read_inputs=True, object1=None, object2=None):
    ''' compute differences in id values between two files, return the missing comoving distance values to check if these are on the shell edge '''
    if read_inputs:
        object1 = read_in_lightcone(input_file1,vars=['id','x','y','z'])
        object2 = read_in_lightcone(input_file2,vars=['id','x','y','z'])
        print('read in lightcone files')
    diff_id1,diff_id2 = compare_ids(object1,object2)
    print('computed differences between id sets')
    mask_id1 = create_mask(object1,diff_id1)
    mask_id2 = create_mask(object2,diff_id2)
    print('created masks')
    missing1 = np.sqrt(object1['x'][mask_id1]**2 + object1['y'][mask_id1]**2 + object1['z'][mask_id1]**2)
    missing2 = np.sqrt(object2['x'][mask_id2]**2 + object2['y'][mask_id2]**2 + object2['z'][mask_id2]**2)
    return missing1,missing2

def histograms(input_file,bins=200, read_inputs=True, object1=None, display_plots=True):
    ''' create projected 2d histograms'''
    if read_inputs:
        object1 = read_in_lightcone(input_file,vars=['x','y','z'])
        print('finished reading in lightcone data')
    plt.figure()
    a1,b1,c1,d1 = plt.hist2d(object1['x'],object1['y'],bins=bins);
    plt.xlabel('x');plt.ylabel('y')
    plt.title(input_file)
    plt.savefig(input_file+'_xy_hist.png',dpi=500)
    plt.figure()
    a2,b2,c2,d2 = plt.hist2d(object1['y'],object1['z'],bins=bins);
    plt.xlabel('y');plt.ylabel('z')
    plt.title(input_file)
    plt.savefig(input_file+'_yz_hist.png',dpi=500)
    plt.figure()
    a3,b3,c3,d3 = plt.hist2d(object1['x'],object1['z'],bins=bins);
    plt.xlabel('x');plt.ylabel('z')
    plt.title(input_file)
    plt.savefig(input_file+'_xz_hist.png',dpi=500)
    if display_plots:
        plt.show()
    return a1,a2,a3

def compare_hists(input_file1,input_file2,bins, read_inputs=True, object1=None, object2=None, display_plots=True):
    a1,a2,a3 = histograms(input_file1,bins=bins, read_inputs=read_inputs, object1=object1, display_plots=display_plots);
    b1,b2,b3 = histograms(input_file2,bins=bins, read_inputs=read_inputs, object1=object2, display_plots=display_plots);
    print('plotting differences')
    plt.figure()
    plt.imshow(a1-b1)
    plt.title('difference plot ')
    plt.xlabel('x');plt.ylabel('y')
    plt.savefig(input_file1+'_'+input_file2+'_xy.png',dpi=500)
    plt.figure()
    plt.imshow(a2-b2)
    plt.title('difference plot')
    plt.xlabel('y');plt.ylabel('z')
    plt.savefig(input_file1+'_'+input_file2+'_yz.png',dpi=500)
    plt.figure()
    plt.imshow(a3-b3)
    plt.title('difference plot')
    plt.xlabel('x');plt.ylabel('z')
    plt.savefig(input_file1+'_'+input_file2+'_xz.png',dpi=500)
    if display_plots:
        plt.show()
    return a1-b1,a2-b2,a3-b3


def compare_vels(input_file1,input_file2,bins=200, read_inputs=True, object1=None, object2=None, display_plots=True):
    ''' compare absolute velocities from two files '''
    if read_inputs:
        object1 = read_in_lightcone(input_file1,vars=['vx','vy','vz'])
        object2 = read_in_lightcone(input_file2,vars=['vx','vy','vz'])
        print('finished reading in lightcone velocities')
    plt.figure()
    plt.hist(np.sqrt(object1['vx']**2+object1['vy']**2+object1['vz']**2),bins=bins,alpha=0.5);
    plt.hist(np.sqrt(object2['vx']**2+object2['vy']**2+object2['vz']**2),bins=bins,alpha=0.5);
    plt.xlabel('v_abs')
    plt.yscale('log')
    plt.savefig('v_'+input_file1+'_'+input_file2+'.png',dpi=500)
    if display_plots:
        plt.show()
    return


def create_healpix_map(input_file, nside=2048,read_inputs=True, object1=None):
    ''' creating healpix map - assuming particle masses are the same. We can remove this assumption easily if we want to check halo maps'''
    map_dens = np.zeros(hp.nside2npix(nside))
    if read_inputs:
        object1 = read_in_lightcone(input_file,vars=['x','y','z'])
        print('read in lightcone values')
    pix_list = hp.vec2pix(nside,object1['x'],object1['y'],object1['z'])
    unique, unique_counts = np.unique(pix_list,return_counts=True)
    map_dens[unique]+= unique_counts
    print('created ngp density map')
    return map_dens


def compare_power(input_file1,input_file2,nside=2048,lmax=2000, read_inputs=True, object1=None, object2 = None, display_plots=True):
    ''' simple power spectrum comparison '''
 
    map_dens1 = create_healpix_map(input_file1,nside,read_inputs=read_inputs, object1=object1)
    map_dens2 = create_healpix_map(input_file2,nside,read_inputs=read_inputs, object1=object2)

    # scaling to remove monopole - this helps with power spectrum computation. 
    mask_octant = compute_mask(nside=nside,octant="upper")
    map_dens1 = scale_map(map_dens1,mask_octant,nside=nside)
    map_dens2 = scale_map(map_dens2,mask_octant,nside=nside)

    print('comparing density maps')
    hp.mollview(map_dens1)
    hp.mollview(map_dens2)
    hp.mollzoom(map_dens1-map_dens2)
    if display_plots:
        plt.show()
    alm_1 = hp.map2alm(map_dens1,lmax=lmax)
    alm_2 = hp.map2alm(map_dens2,lmax=lmax)
    cl_1 = hp.alm2cl(alm_1)
    cl_2 = hp.alm2cl(alm_2)
    print('computed power spectra')
    l = np.arange(lmax+1)
    fsky = 1./8
    plt.figure()
    plt.plot(l,cl_1*l*(l+1)/(2.*np.pi)/fsky)
    plt.plot(l,cl_2*l*(l+1)/(2.*np.pi)/fsky)
    plt.savefig('cls_'+input_file1+'_'+input_file2+'.png',dpi=500)
    if display_plots:
        plt.show()
    plt.figure()
    plt.plot(l,np.abs(cl_1-cl_2)/cl_1)
    plt.savefig('cls_fractional_'+input_file1+'_'+input_file2+'.png',dpi=500)
    if display_plots:
        plt.show()
    return map_dens1,map_dens2,cl_1,cl_2

def check_all(input_file1,input_file2,input_file1_trunc,input_file2_trunc):

    print('reading in lightcone files')
    object1 = read_in_lightcone(input_file1);
    object2 = read_in_lightcone(input_file2);

    ###########
    print('computing comoving distances corresponding to missing id values')
    dist1,dist2 = missing_distances(input_file1,input_file2, read_inputs=False, object1=object1, object2=object2)
    dist1.tofile('dist1_'+input_file1_trunc+'_'+input_file2_trunc+'.txt',sep='\n')
    dist2.tofile('dist2_'+input_file1_trunc+'_'+input_file2_trunc+'.txt',sep='\n')
    print('missing from first file: mean value = ' +str(np.mean(dist1)) +', maximum value = '+str(np.max(dist1))+ ', minimum value = '+str(np.min(dist1))) 
    print(dist1)
    print('missing from second file: mean value = ' +str(np.mean(dist2)) +', maximum value = '+str(np.max(dist2))+ ', minimum value = '+str(np.min(dist2)))
    print(dist2)

    ##########
    print('creating histogram plots')
    diff_xy, diff_yz, diff_xz = compare_hists(input_file1_trunc,input_file2_trunc,bins=np.linspace(0.0,1600.,200), read_inputs=False, object1=object1, object2 = object2, display_plots=False)

    #########
    print('creating velocity plots')
    compare_vels(input_file1_trunc,input_file2_trunc,bins=np.linspace(0.0,10000,200), read_inputs=False, object1=object1, object2=object2, display_plots=False)

    ########
    print('creating healpix maps')
    compare_power(input_file1_trunc,input_file2_trunc,nside=1024,lmax=1000, read_inputs=False, object1=object1, object2 = object2, display_plots=False)
    return;


if __name__ == "__main__":
    input_file1 = sys.argv[1]; input_file2 = sys.argv[2]; input_file1_trunc = sys.argv[3]; input_file2_trunc = sys.argv[4];
    check_all(input_file1,input_file2,input_file1_trunc,input_file2_trunc)
