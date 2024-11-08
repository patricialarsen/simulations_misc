
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
#include "PLHalos.h"
//#include "utils_ngp.h"
#include "pix_funcs.h"

#include "RadiativeCooling.h"

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


// note the difference here is that the halos have a larger extent than the particles so we want them to be assigned to the neighbouring low-res pixels too 
int compute_ranks_count_halos( PLHalos* P, T_Healpix_Base<int> map_lores,  int numranks, vector<int> &send_count, int rank_diff){
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
       fix_arr<int,4> neighbs;
       fix_arr<double,4> weights;
       map_lores.get_interpol(point,  neighbs, weights);
       for (int j=0;j<4;j++){
           int rankn = compute_rank_fullsky(neighbs[j],numranks);
           //rank_list[rankn] = 1;
#ifdef _OPENMP
    omp_set_lock(&(m_lck[rankn]));
#endif
           send_count[rankn] +=1 ;
#ifdef _OPENMP
    omp_unset_lock(&(m_lck[rankn]));
#endif
       }
   }
#ifdef _OPENMP
 for(unsigned i = 0; i < m_lck.size(); i++)omp_destroy_lock(&(m_lck[i]));
#endif
 return 0;
}

void compute_ranks_index_halo( PLHalos* P, T_Healpix_Base<int> map_lores, int numranks, vector<int> send_off, vector<int64_t> &id_particles, int rank_diff){

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
       fix_arr<int,4> neighbs;
       fix_arr<double,4> weights;
       map_lores.get_interpol(point,  neighbs, weights);

       for (int j=0;j<4;j++){
           int ind = neighbs[j];// check order is right   //get_lores_pix(neighbs[j],rank_diff, map_hires);
           int rankn = compute_rank_fullsky(ind,numranks);
           //rank_list[rankn] = 1;
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
#ifdef _OPENMP
 for(unsigned i = 0; i < m_lck.size(); i++)omp_destroy_lock(&(m_lck[i]));
#endif
 return;
}



int assign_sz_cic_halocut(vector<float> &rho, vector<float> &phi, vector<float> &vel,vector<float> &ksz, vector<double> &tsz, PLParticles* P, PLHalos* H,  T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx, float hval, bool borgcube, bool adiabatic, float samplerate, string cloudypath, float masscut, float radiuscut){

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  int64_t npix = map_hires.Npix();
  float pixsize = (4.*3.141529/npix);

  // extra constants 
  const double SIGMAT = 6.65245e-25; // cm^2
  const double MP = 1.672621e-24; // MP in g
  const double MPC_TO_CM = KM_IN_MPC*CM_IN_KM;



  // only used for adiabatic hydro 
  const double MUE = 2.0/(2.0-PRIMORDIAL_Y);
  double NHE = 2.0; // assume helium is fully ionized
  double CHIE = (1.0-PRIMORDIAL_Y*(1.0-NHE/4.0))/(1.0-PRIMORDIAL_Y/2.0);
  const double mu0 = MU_ION;

  // In BorgCube we assumed helium was neutral
  if (borgcube) {
  NHE = 0.0;
  CHIE = (1.0-PRIMORDIAL_Y*(1.0-NHE/4.0))/(1.0-PRIMORDIAL_Y/2.0);
  }

  const double NE_SCALING = (double)CHIE/MUE/MP;
  const double KSZ_CONV = (double)(-SIGMAT*CHIE/MUE/MP/CLIGHT)*(G_IN_MSUN/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize)*hval;
  const double TSZ_CONV = (double)((GAMMA-1.0)*SIGMAT*CHIE/MUE/MELECTRON/CLIGHT/CLIGHT)*(MH/MP)*(G_IN_MSUN/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize)*hval;

   // create halo map 
   // not threaded for now 
   unordered_map<int64_t, float> tag_to_mass;
   unordered_map<int64_t, float> tag_to_radius;
   unordered_map<int64_t, float> tag_to_x;
   unordered_map<int64_t, float> tag_to_y;
   unordered_map<int64_t, float> tag_to_z;


   for (int64_t ii=0; ii<(*H).nparticles; ++ii){
     int64_t halo_tag = (*H).int64_data[0]->at(ii);
     float halo_mass = (*H).float_data[3]->at(ii);
     float halo_radius = (*H).float_data[4]->at(ii);
     float halo_x = (*H).float_data[0]->at(ii);
     float halo_y = (*H).float_data[1]->at(ii);
     float halo_z = (*H).float_data[2]->at(ii);

     tag_to_mass[halo_tag] = halo_mass;
     tag_to_radius[halo_tag] = halo_radius;
     tag_to_x[halo_tag] = halo_x;
     tag_to_y[halo_tag] = halo_y;
     tag_to_z[halo_tag] = halo_z;
   }


  // compute average a value for Radiative Cooling

  double sum_a = 0;
  double sum_a_global;
  int64_t count_a = 0;
  int64_t count_a_global;
  double sum_T =0;
  double sum_mu = 0;
  double sum_mu_global;
  double sum_T_global;
  double min_T = 1.e10;
  double min_mu = 1.e10;
  double min_a = 1.e10;
  double max_T = -1.e10;
  double max_mu = -1.e10;
  double max_a = -1.e10;
  double max_mu_global;
  double max_T_global;
  double max_a_global;
  double min_mu_global;
  double min_T_global;
  double min_a_global;

  // technically this is double work, we need to do this for aa here anyway but could move the rest into the main loop 
  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:sum_a,count_a,sum_mu,sum_T) reduction(max:max_T,max_mu,max_a) reduction(min:min_T,min_mu,min_a)
  #endif
  for (int64_t ii=0; ii<(*P).nparticles; ++ii){
      double aa = (double) (*P).float_data[7]->at(ii);
      double uu = (double) (*P).float_data[10]->at(ii);
      #ifdef HYBRID_SG
      double mu = (double) (*P).float_data[11]->at(ii);
      #else
      double mu = (double) mu0;
      #endif
      uint16_t mask = (*P).mask_data[0]->at(ii);
      double Ti = CONV_U_TO_T*UU_CGS*mu*uu*aa*aa; // UU_CGS is just the km to cm scaling. scale_uu should come from 
      int64_t halo_tag = (*P).int64_data[1]->at(ii);

      if ((isNormGas(mask))&&(tag_to_mass[halo_tag]>masscut)){
      sum_T += Ti;
      sum_mu += mu;
      sum_a += aa;
      min_a = (aa<min_a)?aa:min_a;
      min_T = (Ti<min_T)?Ti:min_T;
      min_mu = (mu<min_mu)?mu:min_mu;
      max_a = (aa>max_a)?aa:max_a;
      max_T = (Ti>max_T)?Ti:max_T;
      max_mu = (mu>max_mu)?mu:max_mu;
      count_a += 1;
      }
  }
  MPI_Reduce(&sum_a,&sum_a_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&sum_T,&sum_T_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&sum_mu,&sum_mu_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&count_a,&count_a_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);

  MPI_Reduce(&min_a,&min_a_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&min_T,&min_T_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&min_mu,&min_mu_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

  MPI_Reduce(&max_a,&max_a_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&max_T,&max_T_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&max_mu,&max_mu_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Bcast(&sum_a_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sum_T_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sum_mu_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&min_a_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&min_T_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&min_mu_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max_a_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max_T_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max_mu_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);


  MPI_Bcast(&count_a_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double aa_av = sum_a_global/(double)count_a_global;
  double T_av = sum_T_global/(double)count_a_global;
  double mu_av = sum_mu_global/(double)count_a_global;

  if (commRank==0){
     cout << "a values: average = "<< aa_av << ", minimum value = "<< min_a_global << ", maximum value = " << max_a_global<< endl;
     cout << "T values (CGS units): average = "<< T_av << ", minimum value = "<< min_T_global << ", maximum value = " << max_T_global<< endl;
     cout << "mu values: average = "<< mu_av << ", minimum value = "<< min_mu_global << ", maximum value = " << max_mu_global<< endl;

  }
  // initialize the cloudy tables

  #ifdef HYBRID_SG
  RadiativeCooling* m_radcool = new RadiativeCooling(cloudypath);
  if (!adiabatic) {
    double Tcmb = 2.725f;
    m_radcool->setTCMB(Tcmb);
    m_radcool->readCloudyScaleFactor((RAD_T)aa_av);
  }
  #endif
  double tsz_tot2 = 0;
  double tsz_tot3 = 0;
  double tsz_tot4 = 0;
  double ne_frac = 0;
  double ne_av = 0;

  int64_t count_part = 0;
  int64_t count_mask = 0;
  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for (int64_t ii=0; ii<(*P).nparticles; ++ii){
      #pragma omp critical
      {
      count_part ++;
      }
      float xd = (float) (*P).float_data[0]->at(ii);
      float yd = (float) (*P).float_data[1]->at(ii);
      float zd = (float) (*P).float_data[2]->at(ii);
      float vx = (float) (*P).float_data[3]->at(ii);
      float vy = (float) (*P).float_data[4]->at(ii);
      float vz = (float) (*P).float_data[5]->at(ii);
      float phid = (float) (*P).float_data[6]->at(ii);
      float aa = (float) (*P).float_data[7]->at(ii);
      double mass = (double) (*P).float_data[9]->at(ii);
      double uu = (double) (*P).float_data[10]->at(ii);
      double mu;
      #ifdef HYBRID_SG
      if (!adiabatic)
        mu = (double) (*P).float_data[11]->at(ii);
      else
        mu = (double) mu0;
      #else
        mu = (double) mu0;
      #endif

      uint16_t mask = (*P).mask_data[0]->at(ii);
      int64_t halo_tag = (*P).int64_data[1]->at(ii);
      if (tag_to_mass[halo_tag]<=masscut){
	     continue; 
      }
      float dist = (xd - tag_to_x[halo_tag])*(xd - tag_to_x[halo_tag]) + (yd - tag_to_y[halo_tag])*(yd - tag_to_y[halo_tag])+(zd - tag_to_z[halo_tag])*(zd - tag_to_z[halo_tag]); 
      dist = dist/(tag_to_radius[halo_tag]*tag_to_radius[halo_tag]);
      if (sqrt(dist)<radiuscut){
             continue; 
      };
      

//      if ((isNormGas(mask))&&(tag_to_mass[halo_tag]>masscut)){

      if (isNormGas(mask)){
        #pragma omp critical
        {
        count_mask++;
        }
      }

      double rhoi = (double) (*P).float_data[12]->at(ii);
      double Vi = mass/rhoi/MPC_IN_CM/MPC_IN_CM/MPC_IN_CM *aa*aa*aa/hval/hval/hval ; //  Vi in physical CGS units
      double Zi = (double) (*P).float_data[13]->at(ii);
      double Yi = (double) (*P).float_data[14]->at(ii);

      int iter = 0;

      rhoi *= G_IN_MSUN*MPC_IN_CM*MPC_IN_CM*MPC_IN_CM/aa/aa/aa*hval*hval;
      RAD_T Xi = (1.0-Yi-Zi);
      RAD_T Ti   = CONV_U_TO_T*UU_CGS*mu*uu*aa*aa; // UU_CGS is just the km to cm scaling 

      //assert(Ti>0);
      /*if ((Ti<=0)&&(isNormGas(mask))){
         cout << "temperature value less than 0 , Ti = "<< Ti << endl;
         cout << "temperature value less than 0 , mu = "<< mu << endl;
         cout << "temperature value less than 0 , uu= "<< uu << endl;

      }*/
      RAD_T nHi  = rhoi*Xi*INV_MH; // density / mass in g
      RAD_T nHIi = 0.0;
      RAD_T nei = 0.0;
      RAD_T mui = (RAD_T)mu;

      if (!adiabatic) {
        #ifdef HYBRID_SG
        if ((Ti>0)&&(isNormGas(mask))){
        RAD_T lambda = (*m_radcool)(Ti, (RAD_T)rhoi, (RAD_T)Zi, (RAD_T)Yi, (RAD_T)aa_av, mui, &iter, false, &nHIi, &nei);
        nei *= nHi;
        }
        #endif
      }
      const double NE_SCALING = (double)CHIE/MUE/MP*G_IN_MSUN/hval;
      double ne_scaled = nei*Vi/(double)NE_SCALING; // this is the SPH volume multiplied by the electron density per unit volume 

      float dist_com2 = xd*xd + yd*yd + zd*zd;
      float vel_los = (vx*xd +vy*yd + vz*zd)/sqrt(dist_com2);

      vec3 vec_val = vec3(xd,yd,zd);
      pointing point = pointing(vec_val);
      fix_arr<int64_t,4> neighbs;
      fix_arr<double,4> weights;
      map_hires.get_interpol(point,  neighbs, weights);
      for (int j=0;j<4;j++){
        int64_t pix_num = neighbs[j];
        int64_t new_idx = -99;
        if (ring_to_idx.count(pix_num)){
          new_idx = ring_to_idx[pix_num];
          assert(new_idx<npix);
          assert(new_idx>=0);
          #pragma omp critical
          {
          if (isNormGas(mask)){//TODO: see below 
          // ideally we should add a flag for SFGas and wind to revert to the adiabatic electron density value
          if (adiabatic){
            tsz[new_idx] += mass*mu*uu/dist_com2*weights[j];   // a^2 factors cancel out in ui and dist_comov2  
            ksz[new_idx] += mass*vel_los/dist_com2/aa*weights[j]; // one factor of a cancels from v_los and dist_comov2         
          }
          else{
            tsz[new_idx] += ne_scaled*mu*uu/dist_com2*weights[j];
            ksz[new_idx] += ne_scaled*vel_los/dist_com2/aa*weights[j]; // one factor of a cancels from v_los and dist_comov2
          }

          tsz_tot2 += ne_scaled*mui*uu/dist_com2*weights[j];
          tsz_tot3 += mass*mu*uu/dist_com2*weights[j];
          tsz_tot4 += mass*mu0*uu/dist_com2*weights[j];
          }
          rho[new_idx] += mass/pixsize*weights[j];
          phi[new_idx] += mass*phid*weights[j];
          vel[new_idx] += mass*vel_los*weights[j]; // mass weighting 
        }
       }
      }
   }
  #ifdef HYBRID_SG
  delete m_radcool;
  #endif
   // for each pixel apply this scaling 
   double tsz_tot = 0;
   for (int64_t j=0; j<ksz.size(); j++){
     tsz[j] = tsz[j]*TSZ_CONV;
     tsz_tot += tsz[j];
     ksz[j] = ksz[j]*KSZ_CONV;
    }
  if (commRank==0){

   cout << "map tsz sum = "<< tsz_tot << endl;
   cout << "particle tsz sum = "<< tsz_tot2*TSZ_CONV << endl;
   cout << "particle tsz sum (adiabatic with read mu ) = "<< tsz_tot3*TSZ_CONV << endl;
   cout << "particle tsz sum (adiabatic with mu0) = "<< tsz_tot4*TSZ_CONV << endl;

}
   return 0;
}

          

void clear_pixel_hydro(int64_t start_idx, int rank_diff, vector<float> &rho , vector<float> &phi, vector<float> &vel, vector<float> &ksz, vector<double> &tsz, vector<double> &xray1, vector<double> &xray2){
    int64_t len_pixel = pow(4,rank_diff);
    for (int64_t idx=start_idx; idx<(len_pixel+start_idx); idx++){
        rho[idx] = 0.0;
        phi[idx] = 0.0;
        vel[idx] = 0.0;
        ksz[idx] = 0.0;
        tsz[idx] = 0.0;
        xray1[idx] = 0.0;
        xray2[idx] = 0.0;
    }
    return;
}

void clear_pixel(int64_t start_idx, int rank_diff, vector<float> &rho , vector<float> &phi, vector<float> &vel){
    int64_t len_pixel = pow(4,rank_diff);
    for (int64_t idx=start_idx; idx<(len_pixel+start_idx); idx++){
        rho[idx] = 0.0;
        phi[idx] = 0.0;
        vel[idx] = 0.0;
    }
    return;
}

