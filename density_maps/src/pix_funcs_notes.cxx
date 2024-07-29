
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
//#include "utils_ngp.h"
#include "pix_funcs.h"

#include "RadiativeCooling.h"

using namespace gio;
using namespace std;


#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t



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




int assign_dm_cic(vector<float> &rho, vector<float> &phi, vector<float> &vel, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx, float hval, float samplerate){

    int64_t npix = map_hires.Npix();
    float pixsize = (4.*3.141529/npix);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
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
	  #pragma omp critical
	  {
          rho[new_idx] += 1./pixsize*weights[j];
          phi[new_idx] += phid*weights[j];
          vel[new_idx] += vel_los*weights[j];
	  }
	}
      }
   }
     return 0;
}



int assign_dm_ngp(vector<float> &rho, vector<float> &phi, vector<float> &vel, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx){

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
        rho[new_idx] += 1./pixsize;
        phi[new_idx] += phid; 
        vel[new_idx] += vel_los;
      }
     }
     return 0;
}

int assign_sz_cic(vector<float> &rho, vector<float> &phi, vector<float> &vel,vector<float> &ksz, vector<float> &tsz, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx, float hval){


  int64_t npix = map_hires.Npix();
  float pixsize = (4.*3.141529/npix);

  const double SIGMAT = 6.65245e-25;
  const double MP = 1.672621e-24; // MP in g
  const double MUE = 2.0/(2.0-PRIMORDIAL_Y); // only for adiabatic, mu should be overwritten otherwise 
  const double NHE    = 0.0; // assume helium is neutral 
  const double CHIE   = (1.0-PRIMORDIAL_Y*(1.0-NHE/4.0))/(1.0-PRIMORDIAL_Y/2.0);
  const double MPC_TO_CM = KM_IN_MPC*CM_IN_KM;
  const double mu0 = MU_ION;
  

  // TODO input the sample rate 
  float samplerate = 1.0; 

  const float NE_SCALING = (float)CHIE/MELECTRON/MP;
  // remember these need to be multiplied by a factor of h later currently 
  const float KSZ_CONV = (float)(-SIGMAT*CHIE/MUE/MP/CLIGHT)*(G_IN_MSUN/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize);
  const float TSZ_CONV = (float)((GAMMA-1.0)*SIGMAT*CHIE/MUE/MELECTRON/CLIGHT/CLIGHT)*(MH/MP)*(G_IN_MSUN/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize);

  string outpath = "/lus/eagle/projects/CosDiscover/nfrontiere/576MPC_RUNS/challenge_problem_576MPC_SEED_1.25e6_NPERH_AGN_2_NEWCOSMO/output/";


  // compute a value 
  double sum_a = 0;
  int64_t count_a = 0;
  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:sum_a,count_a)
  #endif
  for (int64_t ii=0; ii<(*P).nparticles; ++ii){
      double aa = (double) (*P).float_data[7]->at(ii);
      sum_a += aa;
      count_a += 1; 
  }
  MPI_Reduce(&sum_a,&sum_a,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&count_a,&count_a,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&sum_a,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&count_a,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  double aa = sum_a/count_a;

  RadiativeCooling* m_radcool = new RadiativeCooling(outpath);
  double Tcmb = 2.725f;
  m_radcool->setTCMB(Tcmb);
  //double aa = 1.0; // set this properly
  m_radcool->readCloudyScaleFactor((RAD_T)aa);

    /*RAD_T rhoi = rho[nn]*m_hubble*m_hubble*omega_cb*(1.0f/(aa*aa*aa));
    RAD_T mui  = mu[nn];
    RAD_T Ti   = CONV_U_TO_T*UU_CGS*scale_uu*mui*uu[nn]*aa*aa;
    RAD_T Zi   = zmet[nn];
    RAD_T Yi   = yhe[nn];
    RAD_T Xi = (1.0-Yi-Zi);
    RAD_T nHi  = rhoi*RHO_CGS*Xi*INV_MH;
    RAD_T nHIi = 0.0;
    RAD_T nei = 0.0;
    int iter = 0;*/
    //rhoi *= RHO_CGS;
    //RAD_T lambda = (*m_radcool)(Ti, rhoi, Zi, Yi, (RAD_T)aa, mui, &iter, false, &nHIi, &nei);


  //delete m_radcool;


  for (int64_t ii=0; ii<(*P).nparticles; ++ii){

      float xd = (float) (*P).float_data[0]->at(ii);
      float yd = (float) (*P).float_data[1]->at(ii);
      float zd = (float) (*P).float_data[2]->at(ii);
      float vx = (float) (*P).float_data[3]->at(ii);
      float vy = (float) (*P).float_data[4]->at(ii);
      float vz = (float) (*P).float_data[5]->at(ii);
      float phid = (float) (*P).float_data[6]->at(ii);
      float aa = (float) (*P).float_data[7]->at(ii);

      float mass = (float) (*P).float_data[9]->at(ii);
      float uu = (float) (*P).float_data[10]->at(ii);
       
      #ifndef HYBRID_SG
      float mu = (float) (*P).float_data[11]->at(ii);
      #else
      float mu = (float) mu0; 
      #endif 

      float rhoi = (float) (*P).float_data[12]->at(ii);// in PA rho *m_hubble * m_hubble * omega_cb /a^3
      float Zi = (float) (*P).float_data[13]->at(ii); // zi = zmet
      float Yi = (float) (*P).float_data[14]->at(ii);// yi = yhe 
      //if (ii==0){
      int iter = 0;
      aa=1.0;
      //doauble scale_uu = Domain::grid2phys_uu();
 

      //double scale_uu = 1.0;
      //double scale_uu = Domain::grid2phys_uu(); - might need to incorporate unit scalings here

      rhoi *= G_IN_MSUN*MPC_IN_CM*MPC_IN_CM*MPC_IN_CM/aa/aa/aa*hval*hval; 
      RAD_T Xi = (1.0-Yi-Zi);
      //conv_u_to_t converts u to t not including a factor of mu 
      // uu in units of (km/s)^2
      RAD_T Ti   = CONV_U_TO_T*UU_CGS*mu*uu*aa*aa; // UU_CGS is just the km to cm scaling. scale_uu should come from 
      RAD_T nHi  = rhoi*Xi*INV_MH; 
      RAD_T nHIi = 0.0;
      RAD_T nei = 0.0;
      RAD_T lambda = (*m_radcool)(Ti, rhoi, Zi, Yi, (RAD_T)aa, mu, &iter, false, &nHIi, &nei);
      nei *= nHi; // Xi, Yi,Zi all fractions, these should be fine. 
      // mu dimensionless, rho/T we have units. 
      //nei /= (G_IN_MSUN*MPC_IN_CM*MPC_IN_CM*MPC_IN_CM);  // the output nei should be a density in the same units as before 
      //nei currently in physical CGS units

      double Vi = (double)mass/rhoi*G_IN_MSUN/hval; // mass was in Msun/h, needs to be in g, volume is now in (cm)^3
      double ne_scaled = nei*Vi/(double)NE_SCALING; // this is the SPH volume multiplied by the electron density per unit volume 

      //double ne_old = (double)mass;
      // this should be the total number of electrons per particle
      //double ne_old = // chi * rhog / mue /mp - in this case you use rho = mass / vol_pix 
      //nei*mass/rhoi;

      // in this temperature is ( GAMMA = 5/3)
      // uu in comoving non-CGS units, mu dimensionless 
      // (GAMMA-1) MH / KB * mu * uu /a^2  (1.e10 - which is CM_IN_KM^2), assume MH/KB cancel units
 
      // and n is
      // ne = CHIE * rhog / MUE /MP
      // rhog = sum_mass / (A_steradian *  R^2 DR )*a^2
      // ne = CHIE  (sum_i M_i ) /MUE /MP / (AR^2DR)*a^2

      // we then have a prefactor 
      // SIGMAT KB / MELECTRON / C / C (convert to CGS using G_IN_MSUN/ CM_IN_KM^2)


      // and the total gives 
      //
      // DR  (GAMMA-1) MH / KB * mu * uu /a^2 CHIE  (sum_i M_i ) /MUE /MP / (AR^2 DR) *a^2 * SIGMAT KB / MELECTRON / C / C
      // 

      //const float TSZ_CONV = (float)((GAMMA-1.0)*SIGMAT*CHIE/MUE/MELECTRON/CLIGHT/CLIGHT)*(MH/MP)*(G_IN_MSUN/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize);

      // and density 
      // mass * mu * uu / dist_com2   

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
          if (isNormGas(mask)){//TODO: check this
          // ideally we should add a flag for SFGas and wind to revert to the adiabatic electron density value
          ksz[new_idx] += mass*vel_los/dist_com2/aa; // one factor of a cancels from v_los and dist_comov2
          //tsz[new_idx] += mass*mu*uu/dist_com2;   // a^2 factors cancel out in ui and dist_comov2  
          tsz[new_idx] += ne_scaled*mu*uu/dist_com2;
          }
          rho[new_idx] += mass/pixsize*weights[j];
          phi[new_idx] += mass*phid*weights[j];
          vel[new_idx] += mass*vel_los*weights[j]; // mass weighting 
        }
      }
   }
  delete m_radcool;
   // for each pixel apply this scaling 
   for (int64_t j=0; j<ksz.size(); j++){
     tsz[j] = tsz[j]*TSZ_CONV;
     ksz[j] = ksz[j]*KSZ_CONV;
    }
   return 0;
}




int assign_sz_ngp(vector<float> &rho, vector<float> &phi, vector<float> &vel,vector<float> &ksz, vector<float> &tsz, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx){

    int64_t npix = map_hires.Npix();
    float pixsize = (4.*3.141529/npix);

  const double SIGMAT = 6.65245e-25;
  const double MP = 1.672621e-24; // MP in g
  const double MUE = 2.0/(2.0-PRIMORDIAL_Y); // only for adiabatic, mu should be overwritten otherwise 
  const double NHE    = 0.0; // assume helium is neutral 
  const double CHIE   = (1.0-PRIMORDIAL_Y*(1.0-NHE/4.0))/(1.0-PRIMORDIAL_Y/2.0);
  const double MPC_TO_CM = KM_IN_MPC*CM_IN_KM;
  const double mu0 = MU_ION;

  // TODO input the sample rate 
     float samplerate = 1.0;
  //
  const float KSZ_CONV = (float)(-SIGMAT*CHIE/MUE/MP/CLIGHT)*(G_IN_MSUN/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize);
  const float TSZ_CONV = (float)((GAMMA-1.0)*SIGMAT*CHIE/MUE/MELECTRON/CLIGHT/CLIGHT)*(MH/MP)*(G_IN_MSUN/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize);



    for (int64_t ii=0; ii<(*P).nparticles; ++ii){
      float xd = (float) (*P).float_data[0]->at(ii);
      float yd = (float) (*P).float_data[1]->at(ii);
      float zd = (float) (*P).float_data[2]->at(ii);
      float vx = (float) (*P).float_data[3]->at(ii);
      float vy = (float) (*P).float_data[4]->at(ii);
      float vz = (float) (*P).float_data[5]->at(ii);
      float phid = (float) (*P).float_data[6]->at(ii);
      // in the case of hydro 
      float aa = (float) (*P).float_data[7]->at(ii);

      float mass = (float) (*P).float_data[9]->at(ii);
      float uu = (float) (*P).float_data[10]->at(ii);
      #ifndef HYBRID_SG
      float mu = (float) (*P).float_data[11]->at(ii);
      #else
      float mu = (float) mu0; 
      #endif


      // hydro parameters - we are missing some here - rho, yhe, zmet
      // we might want to use rho for another density measurement 

      float dist_com2 = xd*xd + yd*yd + zd*zd;
      float vel_los = (vx*xd +vy*yd + vz*zd)/sqrt(dist_com2);

      vec3 vec_val = vec3(xd,yd,zd);
      int64_t pix_num = map_hires.vec2pix(vec_val);
      int64_t new_idx = -99;
      if (ring_to_idx.count(pix_num)){
        new_idx = ring_to_idx[pix_num];
        if (isNormGas(mask)){ // double check this 
        ksz[new_idx] += mass*vel_los/dist_com2/aa; // one factor of a cancels from v_los and dist_comov2
        tsz[new_idx] += mass*mu*uu/dist_com2;   // a^2 factors cancel out in ui and dist_comov2  
        }
        rho[new_idx] += 1./pixsize;
        phi[new_idx] += phid;
        vel[new_idx] += vel_los;
      }
     }
     for (int64_t j=0; j<(npix*pixnum_start.size()); j++){
       tsz[j] = tsz[j]*TSZ_CONV;
       ksz[j] = ksz[j]*KSZ_CONV;
     }
 
     return 0;
}

void initialize_pixel_hydro(int pix_val,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, vector<float> &rho , vector<float> &phi, vector<float> &vel, vector<float> &ksz, vector<float> &tsz, int64_t &count, vector<int64_t> &start_idx, vector<int64_t> &end_idx, vector<int64_t> &pixnum_start, vector<int64_t> &pixnum_end, int rank_diff,  unordered_map<int64_t, int64_t>* ring_to_idx){
    int64_t npix = map_hires.Npix();
    
    start_idx.push_back(count);
    
    int64_t start_count = count;
    int64_t pixx = (int64_t) pix_val;
    
    int64_t start_pixel, len_pixel; 
    len_pixel = pow(4,rank_diff);
    start_pixel = pix_val*len_pixel;
    pixnum_start.push_back(start_pixel);
    pixnum_end.push_back(start_pixel+len_pixel);
    for (int64_t idx=0; idx<len_pixel; idx++){
                        int64_t pix_ring = map_hires.nest2ring(start_pixel+idx);
                        (*ring_to_idx)[pix_ring] = count;
                        rho.push_back(0.0); // maybe initialize these to full size first based on number of pixels
                        phi.push_back(0.0);
                        vel.push_back(0.0);
                        ksz.push_back(0.0);
                        tsz.push_back(0.0);
                         count++;
    }
     end_idx.push_back(count);

    return;
}



void initialize_pixel(int pix_val,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, vector<float> &rho , vector<float> &phi, vector<float> &vel, int64_t &count, vector<int64_t> &start_idx, vector<int64_t> &end_idx, vector<int64_t> &pixnum_start, vector<int64_t> &pixnum_end, int rank_diff,  unordered_map<int64_t, int64_t>* ring_to_idx){
    // initialize large-scale pixel to zero values
    int64_t npix = map_hires.Npix();

    start_idx.push_back(count);

    int64_t start_count = count;
    int64_t pixx = (int64_t) pix_val;

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
