
#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#if defined(_OPENMP) 
#include <omp.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <stdint.h>

#include "chealpix.h"
#include "GenericIO.h"
#include "BasicDefinition.h"


#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>

#include <vector>

#include <climits>
#include <algorithm>

#include "healpix_base.h"


#include "pointing.h"
#include "vec3.h"

#include "PLParticles.h"
#include "pix_funcs_gravonly.h"


using namespace gio;
using namespace std;


#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t
//#define RAD_T float


// compute_rank_fullsky - round robin assignation from full sky assumption 
// compute_rank_partsky - round robin assignation from pixel list 
// get_lores_pix - get low resolution pixel associated with high res one 
// compute_ranks_count - compute number of pixels to send to each rank for communication counts
// compute_ranks_index - compute index of pixels to send to each rank
// get_pix_list_rank - get list of low res pixels in rank using above functions
// initialize_pixel - initialize vector of map on rank to 0 values 
// assign_dm_ngp - use NGP interpolation to add values to map 
// assign_dm_cic - use CIC interpolation to add values to map 


int compute_rank_fullsky(int pix_val, int numranks){
    // assign particles based on pixel value 'pix_val' a rank value from 'numranks' mpi ranks 
    	int rankn = pix_val%numranks;
	return rankn;
}


int compute_rank_partsky(int pix_val, int numranks, vector<int> pix_list){
   // assign particles based on pixel value 'pix_val' a rank value from 'numranks' mpi ranks, given a list of pixel values 'pix_list'
   int rankn=-99;
   int i;
   for (i=0;i<pix_list.size();++i){
       if (pix_list[i]==pix_val){
       rankn = i%numranks;
       return rankn;
      }
   }
  assert(rankn!=-99);
  return -1;
}

int get_lores_pix(int64_t pix_val, int rank_diff, T_Healpix_Base<int64_t> map_hires){
     int64_t pix_val_nest = map_hires.ring2nest(pix_val); // if in ring order 
     int64_t len_pixel = pow(4,rank_diff);
     int npix_low = (int)((float)pix_val_nest/(float)len_pixel); // you want this to truncate
     return npix_low; 
}




int compute_ranks_count( PLParticles* P, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int numranks, vector<int> &send_count, int rank_diff){
#ifdef _OPENMP
 std::vector<omp_lock_t> m_lck;
 m_lck.resize(numranks);
 for(unsigned i = 0; i < m_lck.size(); i++)omp_init_lock(&(m_lck[i]));
#pragma omp parallel for
#endif
     for (int64_t k=0;k<(*P).nparticles;k++){
       float xx =  (*P).float_data[0]->at(k); 
       float yy =  (*P).float_data[1]->at(k);
       float zz =  (*P).float_data[2]->at(k); // assign to x,y,z listed naems 

       vec3 vec_val = vec3(xx,yy,zz);
       pointing point = pointing(vec_val);
       fix_arr<int64_t,4> neighbs;
       fix_arr<double,4> weights;
       map_hires.get_interpol(point,  neighbs, weights); 
       unordered_map<int,int> rank_list;
       for (int j=0;j<4;j++){
           int ind = get_lores_pix(neighbs[j],rank_diff,map_hires);
	   // change this to be based on a list 
	   int rankn = compute_rank_fullsky(ind,numranks);//,pix_list);
	   if (rank_list.count(rankn)){
		   // do nothing 
           }
	   else{
           rank_list[rankn] = 1;
#ifdef _OPENMP
    omp_set_lock(&(m_lck[rankn]));
#endif
           send_count[rankn] +=1 ; 
#ifdef _OPENMP
    omp_unset_lock(&(m_lck[rankn]));
#endif

	   }
       }
   }
#ifdef _OPENMP
 for(unsigned i = 0; i < m_lck.size(); i++)omp_destroy_lock(&(m_lck[i]));
#endif
 return 0;
}



void compute_ranks_index( PLParticles* P, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int numranks, vector<int> send_off, vector<int64_t> &id_particles, int rank_diff){

   vector<int64_t> count_rank(numranks,0);
#ifdef _OPENMP
 std::vector<omp_lock_t> m_lck;
 m_lck.resize(numranks);
 for(unsigned i = 0; i < m_lck.size(); i++)omp_init_lock(&(m_lck[i]));
#pragma omp parallel for
#endif
   for (int64_t k=0;k<(*P).nparticles;k++){
       float xx =  (*P).float_data[0]->at(k);
       float yy =  (*P).float_data[1]->at(k);
       float zz =  (*P).float_data[2]->at(k); // assign to x,y,z listed names 

       vec3 vec_val = vec3(xx,yy,zz);
       pointing point = pointing(vec_val);
       fix_arr<int64_t,4> neighbs;
       fix_arr<double,4> weights;
       map_hires.get_interpol(point,  neighbs, weights);
       unordered_map<int, int> rank_list;

       for (int j=0;j<4;j++){
           int ind = get_lores_pix(neighbs[j],rank_diff, map_hires);
           int rankn = compute_rank_fullsky(ind,numranks);
	   if (rank_list.count(rankn)){
	   }
	   else{
           rank_list[rankn] = 1;
#ifdef _OPENMP
    omp_set_lock(&(m_lck[rankn]));
#endif
           id_particles[send_off[rankn]+count_rank[rankn]] = k;
	   count_rank[rankn]++;
#ifdef _OPENMP
    omp_unset_lock(&(m_lck[rankn]));
#endif
	   }
       }
   }
#ifdef _OPENMP
 for(unsigned i = 0; i < m_lck.size(); i++)omp_destroy_lock(&(m_lck[i]));
#endif
 return;
}




void get_pix_list_rank(int octant, int rank, int numranks, int64_t npix_lores, vector<int> pix_list_oct, vector<int> pix_list_noct,  vector<int> &lores_pix){
    // get list of low resolution pixels in the rank
    // octant = 1 is octant list (need boundary pixels defined) 
    // octant = 0 is full sky list 
    // octant = 2 is given list 

    if (octant==1){
        for (int i=0;i<pix_list_oct.size()+pix_list_noct.size();++i){
            if (i%numranks==rank){
                if (i<pix_list_oct.size()){
                    // pixel is in octant
                    lores_pix.push_back(pix_list_oct[i]);
                }
                else{
                    // include pixels not in octant
                    lores_pix.push_back(pix_list_noct[i-pix_list_oct.size()]);
                }
            }
        }
    }
    else if (octant==2){
	// given pixel list, round robin 
        for (int i=0;i<pix_list_oct.size();++i){
            if (i%numranks==rank){
                    lores_pix.push_back(pix_list_oct[i]);
                }
 	}
    }
    else{
        for (int i=0;i<npix_lores;++i){
            if (i%numranks==rank){
                lores_pix.push_back(i);
            }
        }
    }
    return;
}

int assign_dm_cic(vector<double> &rho, vector<double> &phi, vector<double> &vel, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx, float samplerate){

    int64_t npix = map_hires.Npix();
    float pixsize = (4.*3.141529/npix);
    for (int64_t ii=0; ii<(*P).nparticles; ++ii){
      float xd = (float) (*P).float_data[0]->at(ii);
      float yd = (float) (*P).float_data[1]->at(ii);
      float zd = (float) (*P).float_data[2]->at(ii);
      float vx = (float) (*P).float_data[3]->at(ii);
      float vy = (float) (*P).float_data[4]->at(ii);
      float vz = (float) (*P).float_data[5]->at(ii);
      float phid = (float) (*P).float_data[6]->at(ii);

      float dist_com2 = xd*xd + yd*yd + zd*zd;
      float vel_los = (vx*xd +vy*yd + vz*zd)/sqrt(dist_com2);

      vec3 vec_val = vec3(xd,yd,zd);
      pointing point = pointing(vec_val);
      fix_arr<int64_t,4> neighbs;
      fix_arr<double,4> weights;
      map_hires.get_interpol(point,  neighbs, weights);

      // each time a particle has a neigbour on the rank we add the contribution
      for (int j=0;j<4;j++){
        int64_t pix_num = neighbs[j];
        int64_t new_idx = -99;
        if (ring_to_idx.count(pix_num)){
	  new_idx = ring_to_idx[pix_num];
          rho[new_idx] += (double) (1./pixsize*weights[j]/samplerate);
          phi[new_idx] += (double) (phid*weights[j]/samplerate);
          vel[new_idx] += (double) (vel_los*weights[j]/samplerate);
	}
      }
   }
     return 0;
}



int assign_dm_ngp(vector<double> &rho, vector<double> &phi, vector<double> &vel, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx){

    int64_t npix = map_hires.Npix();
    float pixsize = (4.*3.141529/npix);
    for (int64_t ii=0; ii<(*P).nparticles; ++ii){
      float xd = (float) (*P).float_data[0]->at(ii);
      float yd = (float) (*P).float_data[1]->at(ii);
      float zd = (float) (*P).float_data[2]->at(ii);
      float vx = (float) (*P).float_data[3]->at(ii);
      float vy = (float) (*P).float_data[4]->at(ii);
      float vz = (float) (*P).float_data[5]->at(ii);
      float phid = (float) (*P).float_data[6]->at(ii);
      
      float dist_com2 = xd*xd + yd*yd + zd*zd;
      float vel_los = (vx*xd +vy*yd + vz*zd)/sqrt(dist_com2);

      vec3 vec_val = vec3(xd,yd,zd);
      int64_t pix_num = map_hires.vec2pix(vec_val);
      int64_t new_idx = -99;
      if (ring_to_idx.count(pix_num)){
        new_idx = ring_to_idx[pix_num];
        rho[new_idx] += (double)(1./pixsize);
        phi[new_idx] += (double) phid; 
        vel[new_idx] += (double) vel_los;
      }
     }
     return 0;
}




void clear_pixel(int64_t start_idx, int rank_diff, vector<double> &rho , vector<double> &phi, vector<double> &vel){
    int64_t len_pixel = pow(4,rank_diff);
    for (int64_t idx=start_idx; idx<(len_pixel+start_idx); idx++){
        rho[idx] = 0.0;
        phi[idx] = 0.0;
        vel[idx] = 0.0;
    }
    return;
}




void initialize_pixel(int pix_val,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, vector<double> &rho , vector<double> &phi, vector<double> &vel, int64_t &count, vector<int64_t> &start_idx, vector<int64_t> &end_idx, vector<int64_t> &pixnum_start, vector<int64_t> &pixnum_end, int rank_diff,  unordered_map<int64_t, int64_t>* ring_to_idx){
    // initialize large-scale pixel to zero values
    int64_t npix = map_hires.Npix();

    start_idx.push_back(count);

    int64_t start_count = count;

    int64_t start_pixel, len_pixel; 
    len_pixel = pow(4,rank_diff);
    start_pixel = pix_val*len_pixel;
    pixnum_start.push_back(start_pixel);
    pixnum_end.push_back(start_pixel+len_pixel);
    for (int64_t idx=0; idx<len_pixel; idx++){
	                int64_t pix_ring = map_hires.nest2ring(start_pixel+idx);
			(*ring_to_idx)[pix_ring] = count;
                        rho.push_back(0.0);
                        phi.push_back(0.0);
                        vel.push_back(0.0); 
                         count++;
    }
     end_idx.push_back(count);

    return;
}
