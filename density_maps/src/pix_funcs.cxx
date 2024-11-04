
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

//#include "RadiativeCooling.h"

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
     int64_t npix_low_tmp = pix_val_nest/len_pixel; // you want this to truncate
     int npix_low = (int)npix_low_tmp;
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



/*
int check_xray_halo( PLParticles* P, float hval, bool borgcube, bool adiabatic,  float samplerate, string cloudypath){

  int commRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

  // extra constants 
  const double SIGMAT = 6.65245e-25; // cm^2
  const double MP = 1.672621e-24; // MP in g
  const double MPC_TO_CM = KM_IN_MPC*CM_IN_KM;
  // only used for adiabatic hydro 
  const double MUE = 2.0/(2.0-PRIMORDIAL_Y);
  double NHE = 2.0; // assume helium is fully ionized
  double CHIE = (1.0-PRIMORDIAL_Y*(1.0-NHE/4.0))/(1.0-PRIMORDIAL_Y/2.0);
  const double mu0 = MU_ION;

  const double XRAY_CONV = (double) (1.0/samplerate)*(1.0/pixsize)/(4.0*3.141529);
  double aa_av = 1.0;
  #ifdef HYBRID_SG
  RadiativeCooling* m_radcool = new RadiativeCooling(cloudypath);
  if (!adiabatic) {
    double Tcmb = 2.725f;
    m_radcool->setTCMB(Tcmb);
    m_radcool->readCloudyScaleFactor((RAD_T)aa_av);
    if (commRank==0){
    cout << "Initialized cloudy tables"<<endl;
    }


  }
  #endif

  int64_t count_part = 0;
  int64_t count_mask = 0;
  double xray1=0;
  double xray2=0;

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
      float aa = 1.0; 
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

      double y = Yi/Xi; // nHe/nH * (mHe/mH)
      double Yp = y/(1.0+y); // Metal-free helium fraction

      RAD_T nHi  = rhoi*Xi*INV_MH; // density / mass in g
      RAD_T nHIi = 0.0;
      RAD_T nei = 0.0;
      RAD_T mui = (RAD_T)mu;
      RAD_T LCval = 0.0;
      RAD_T LHval = 0.0;
      RAD_T mue = 0.0;
      double Li_free_Bolo = 0.0;
      double Li_free_ROSAT = 0.0;

      if (!adiabatic) {
        #ifdef HYBRID_SG
        if ((Ti>0)&&(isNormGas(mask))){
        // change to this
        RAD_T lambda = (*m_radcool)(Ti, (RAD_T)rhoi, (RAD_T)Zi, (RAD_T)Yp, (RAD_T)aa_av, mui, &iter, false, &nHIi, &nei, &LCval, &LHval, &mue);

        //RAD_T lambda = (*m_radcool)(Ti, (RAD_T)rhoi, (RAD_T)Zi, (RAD_T)Yi, (RAD_T)aa_av, mui, &iter, false, &nHIi, &nei);
        nei *= nHi;
        }

        // PL: I'm going to assume this mask check is fine 
        // we're using aa here, in the SZ stuff we used aa_av as the average value as that's where we're setting the cloudy redshift to be
        //

        if(XrayGasCondition(mask, Ti)){
          //Calculate L500:
          double mCGS = mass*G_IN_MSUN/hval;//Mass in grams
          double Li_free = 0.0;

          // I've confirmed rhoi is equivalent to rhoCGS
          PRIMORDIAL_FREE_FREE_EMISSIVITY(Li_free, Ti, rhoi, mCGS);

          double Li_free_analytic = 0.0;
          double Xe = 1.16;
          double Xi = 1.079;
          double Z = sqrt(1.15);
          double gff = 1.3;
          double planck = 6.6261e-27; // in cm^2*g/s

          double prefactor = 6.8e-38*KB_CONST/planck*Z*Z*gff*Xe*Xi/(Xi+Xe)/(Xi+Xe); // units Erg/K cm^-2 g^-1 s 
          // units here are confusing
          double scaling = sqrt(Ti)*(rhoi/mui/MH)*(rhoi/mui/MH)*Vi;  // units K^1/2 cm^-3
                                                                     // MH rather than MP because of mui definition
          Li_free_analytic = prefactor * scaling;
          // they then claim a erg K^1/2 cm^3 in Pakmor which I guess gives 
          // erg^2  s cm^-2/g units in total. 
          //
          // luminosity should be in erg/s
          // an erg is g cm^2/s^2
          // so erg  * g cm^2 /s^2 * s *cm^-2 /g
          // gives erg/s 

          //Apply frequency band 
          double Bolo_int = exp(-XRAY_BAND_BOLO_EMIN/(KB_KEV*Ti*aa)) - exp(-XRAY_BAND_BOLO_EMAX/(KB_KEV*Ti*aa));
          //Integrate bolometric band (note the 1/aa scalefactor for redshift dependence)
          double ROSAT_int = exp(-XRAY_BAND_ROSAT_EMIN/(KB_KEV*Ti*aa)) - exp(-XRAY_BAND_ROSAT_EMAX/(KB_KEV*Ti*aa));
          //Integrate ROSAT band (note the 1/aa scalefactor for redshift dependence)

          Li_free_Bolo = Li_free*Bolo_int;//Li_free*Bolo_int;//Bolometric Luminosity
          Li_free_ROSAT = Li_free*ROSAT_int;//ROSAT Luminosity
                                                   // for now let's swap the Bolo luminosity for the analytic one
        }
        #endif
      }

          #pragma omp critical
          {
          if (isNormGas(mask)){//TODO: see below 
            xray1 += Li_free_Bolo;///dist_com2*MPC_IN_CM*MPC_IN_CM*aa*aa*hval*hval;
            xray2 += Li_free_ROSAT;///dist_com2*MPC_IN_CM*MPC_IN_CM*aa*aa*hval*hval; // factor of a^2 from luminosity distance             
        }
       }
      }
   //}

  if (commRank==0){
  cout << "particle count is = "<< count_part << endl;
  cout << "masked particle count is = "<< count_mask << endl;
  cout << "xray Bolometric total is = "<< xray1 << endl;
  cout << "xray ROSAT total is = "<< xray2<<endl;
  }

  #ifdef HYBRID_SG
  delete m_radcool;
  #endif

   return 0;


}
*/



int assign_sz_xray_cic(vector<double> &rho, vector<double> &phi, vector<double> &vel,vector<double> &ksz, vector<double> &tsz, vector<double> &xray1, vector<double> &xray2, vector<double> &xray3, vector<double> &xray4, vector<double> &temp, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx, float hval, bool borgcube, bool adiabatic, float samplerate, string cloudypath){

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
  const double XRAY_CONV = (double) (1.0/samplerate)*(1.0/pixsize)/(4.0*3.141529);

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

      if (isNormGas(mask)){
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
  RadiativeCooling* m_radcool = new RadiativeCooling(cloudypath);
  m_radcool->setTCMB(2.2725f);
  m_radcool->readCloudyTable(min_a_global,max_a_global);

    if (commRank==0){
    cout << "Initialized cloudy tables (new version) at "<< cloudypath <<endl;
    }

  RAD_T m_al, m_ah;
  RAD_T m_Ypmin, m_dYp;
  RAD_T m_Rmin, m_dR;
  RAD_T m_nHmin, m_dnH;
  RAD_T m_Tmin, m_dT;
  int m_nYp;
  int m_nR;
  int m_nnH;
  int m_nT;
  int64_t m_tsize;
  vector<RAD_T>* lcTable_l;
  vector<RAD_T>* lhTable_l;
  vector<RAD_T>* feTable_l;
  vector<RAD_T>* fnTable_l;
  vector<RAD_T>* lcTable_h;
  vector<RAD_T>* lhTable_h;
  vector<RAD_T>* feTable_h;
  vector<RAD_T>* fnTable_h;
  vector<RAD_T>* X0Table_l;
  vector<RAD_T>* X1Table_l;
  vector<RAD_T>* X2Table_l;
  vector<RAD_T>* X3Table_l;
  vector<RAD_T>* X0Table_h;
  vector<RAD_T>* X1Table_h;
  vector<RAD_T>* X2Table_h;
  vector<RAD_T>* X3Table_h;

  // TODO: check check check check. I think this should come in from the header but we'll see
  RAD_T* X0TablePtr_l;
  RAD_T* X0TablePtr_h;
  RAD_T* X1TablePtr_l;
  RAD_T* X1TablePtr_h;
  RAD_T* X2TablePtr_l;
  RAD_T* X2TablePtr_h;
  RAD_T* X3TablePtr_l;
  RAD_T* X3TablePtr_h;

  // note: need to set m_al an m_ah correctly for scale factor interpolation
  m_radcool->getTableParameters(m_al, m_ah, m_Ypmin, m_dYp, m_Rmin, m_dR, m_nHmin, m_dnH, m_Tmin, m_dT,
                         m_nYp, m_nR, m_nnH, m_nT, m_tsize);
  m_radcool->getTableData(&lcTable_l, &lhTable_l, &feTable_l, &fnTable_l, &lcTable_h, &lhTable_h, &feTable_h, &fnTable_h,
                               &X0Table_l, &X1Table_l, &X2Table_l, &X3Table_l, &X0Table_h, &X1Table_h, &X2Table_h, &X3Table_h);
  int64_t tabsize = int64_t(m_nYp)*int64_t(m_nR)*int64_t(m_nnH)*int64_t(m_nT);
  assert(tabsize == m_tsize);//sanity check

  if(m_radcool != NULL) { //TODO: I think this mask is not needed, we are not at the individual particle level yet
    X0TablePtr_l = X0Table_l->data(); // these tables are now in the header
    X0TablePtr_h = X0Table_h->data();
    X1TablePtr_l = X1Table_l->data();
    X1TablePtr_h = X1Table_h->data();
    X2TablePtr_l = X2Table_l->data();
    X2TablePtr_h = X2Table_h->data();
    X3TablePtr_l = X3Table_l->data();
    X3TablePtr_h = X3Table_h->data();

  }


  //if (!adiabatic) {
    //double Tcmb = 2.725f;
    //m_radcool->setTCMB(Tcmb); // no longer needed?
    //m_radcool->readCloudyScaleFactor((RAD_T)aa_av); // TODO: now a table read, can use min_a_global, max_a_global
    if (commRank==0){
    cout << "Initialized cloudy tables (new version) "<<endl;
    }
    

  //}
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

      double y = Yi/Xi; // nHe/nH * (mHe/mH)
      double Yp = y/(1.0+y); // Metal-free helium fraction
 
      // UU_CGS is CM_IN_KM**2
     // U_TO_T is ~ MH/KB units of g* K/Erg  (uu units?)
     // let's assume T in K here 

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
      RAD_T LCval = 0.0;
      RAD_T LHval = 0.0;
      RAD_T mue = 0.0;
      double Li_free_Bolo = 0.0;
      double Li_free_ROSAT = 0.0; 
      double Li_free_ErositaLo = 0.0;
      double Li_free_ErositaHi = 0.0;


      if (!adiabatic) {
        #ifdef HYBRID_SG
        if ((Ti>0)&&(isNormGas(mask))){
	// change to this
	
	// TODO: Q - do I need to update the CLOUDY calls for SZ as well?
	RAD_T lambda = (*m_radcool)(Ti, (RAD_T)rhoi, (RAD_T)Zi, (RAD_T)Yp, (RAD_T)aa_av, mui, &iter, false, &nHIi, &nei, &LCval, &LHval, &mue);

        //TODO: possibly new - don't think so?
         //RAD_T lambda = -(*m_radcool)(RAD_T(Ti), RAD_T(rhoInput), RAD_T(zmetp), RAD_T(Yp), RAD_T(aa), mup, &iter, false, &nHIi, &nei, &LCval, &LHval, &mue);//call returns nHIi and ne in 1/nH units
#ifdef UVCONVERT
//If UVCONVERT IS ON Lambda_code = lambda_CGS * KM_IN_MPC * UU_SI/RHO_CGS... So we want to make lambda in CGS unnits and apply conversions below
         double convFactor = RHO_CGS/(KM_IN_MPC*UU_SI);
         lambda *= convFactor;
         LCval *= convFactor;
         LHval *= convFactor;
#endif
//TODO: new stuff 
//#ifdef UV_CLOUDY_TALBES
	     double Ri = X_SOLAR_OVER_Z_SOLAR * (Zi/Xi);
             //double Ri = X_SOLAR_OVER_Z_SOLAR * (zmetp/Xi); // zmetp is Zi in ours 
             float Yp0 = (Yp-m_Ypmin)/m_dYp;
             float R0 = (Ri-m_Rmin)/m_dR;
             float nH0 = (log10(nHi)-m_nHmin)/m_dnH;
             float T0 = (log10(Ti)-m_Tmin)/m_dT;
//#endif

        //RAD_T lambda = (*m_radcool)(Ti, (RAD_T)rhoi, (RAD_T)Zi, (RAD_T)Yi, (RAD_T)aa_av, mui, &iter, false, &nHIi, &nei);
        nei *= nHi;
	}
	
	// PL: I'm going to assume this mask check is fine 
	// we're using aa here, in the SZ stuff we used aa_av as the average value as that's where we're setting the cloudy redshift to be
	//
	
        if(XrayGasCondition(mask, Ti)){
          //Calculate L500:
          double mCGS = mass*G_IN_MSUN/hval;//Mass in grams
          //double Li_free_Bolo = 0.0;
	  //double Li_free_ROSAT = 0.0;

	  double Ri = X_SOLAR_OVER_Z_SOLAR * (Zi/Xi);
          //double Ri = X_SOLAR_OVER_Z_SOLAR * (zmetp/Xi); // zmetp is Zi in ours
          float Yp0 = (Yp-m_Ypmin)/m_dYp;
          float R0 = (Ri-m_Rmin)/m_dR;
          float nH0 = (log10(nHi)-m_nHmin)/m_dnH;
          float T0 = (log10(Ti)-m_Tmin)/m_dT;

	   
	  // I've confirmed rhoi is equivalent to rhoCGS
          CLOUDY_TAB5D(Li_free_Bolo, X3TablePtr_l, X3TablePtr_h, m_al, m_ah, aa, Yp0, R0, nH0, T0, m_nYp, m_nT, m_nnH, m_nR);
          Li_free_Bolo = Vi*nHi*nHi*pow(10.0, Li_free_Bolo); // Table is in log10(L/nH^2) units
          CLOUDY_TAB5D(Li_free_ROSAT, X0TablePtr_l, X0TablePtr_h, m_al, m_ah, aa, Yp0, R0, nH0, T0, m_nYp, m_nT, m_nnH, m_nR);
          Li_free_ROSAT = Vi*nHi*nHi*pow(10.0, Li_free_ROSAT); // Table is in log10(L/nH^2) units
          CLOUDY_TAB5D(Li_free_ErositaLo, X1TablePtr_l, X1TablePtr_h, m_al, m_ah, aa, Yp0, R0, nH0, T0, m_nYp, m_nT, m_nnH, m_nR);
          Li_free_ErositaLo = Vi*nHi*nHi*pow(10.0, Li_free_ErositaLo); // Table is in log10(L/nH^2) units
          CLOUDY_TAB5D(Li_free_ErositaHi, X2TablePtr_l, X2TablePtr_h, m_al, m_ah, aa, Yp0, R0, nH0, T0, m_nYp, m_nT, m_nnH, m_nR);
          Li_free_ErositaHi = Vi*nHi*nHi*pow(10.0, Li_free_ErositaHi); // Table is in log10(L/nH^2) units


        }
        #endif
      } 

      const double NE_SCALING = (double)CHIE/MUE/MP*G_IN_MSUN/hval; 
      double ne_scaled = nei*Vi/(double)NE_SCALING; // this is the SPH volume multiplied by the electron density per unit volume 
      // double ne_scaled = mass/(mue*MH); // check that mH is in the right units. probably an h factor difference here 
      //double ne_scaled = mass/(mue*MH)/CHIE*MUE*MP;//      //ne*Vi is in units of nH (per some unit volume) * volume
						 //      mass /mue*MH is in units of total number. Assuming units match

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
	    xray4[new_idx] += Li_free_Bolo/dist_com2*MPC_IN_CM*MPC_IN_CM*aa*aa*hval*hval*weights[j];
            xray1[new_idx] += Li_free_ROSAT/dist_com2*MPC_IN_CM*MPC_IN_CM*aa*aa*hval*hval*weights[j]; // factor of a^2 from luminosity distance		    
            xray2[new_idx] += Li_free_ErositaLo/dist_com2*MPC_IN_CM*MPC_IN_CM*aa*aa*hval*hval*weights[j];
            xray3[new_idx] += Li_free_ErositaHi/dist_com2*MPC_IN_CM*MPC_IN_CM*aa*aa*hval*hval*weights[j]; // factor of a^2 from luminosity distance        
	    temp[new_idx]  += Ti*Li_free_Bolo/dist_com2*MPC_IN_CM*MPC_IN_CM*aa*aa*hval*hval*weights[j];  
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
  if (commRank==0){

  cout << "particle count is = "<< count_part << endl;
  cout << "masked particle count is = "<< count_mask << endl;
  }
  #ifdef HYBRID_SG
  delete m_radcool;
  #endif
   // for each pixel apply this scaling 
   int64_t tsz_nzero = 0;
   int64_t ksz_nzero = 0;
   int64_t xray1_nzero = 0;
   int64_t xray2_nzero = 0;
   int64_t xray3_nzero = 0;
   int64_t xray4_nzero = 0;
   int64_t temp_nzero = 0;
   int64_t rho_nzero = 0;
   int64_t phi_nzero = 0;
   int64_t vel_nzero = 0;

   int64_t tsz_nzero_global = 0;
   int64_t ksz_nzero_global = 0;
   int64_t xray1_nzero_global = 0;
   int64_t xray2_nzero_global = 0;
   int64_t xray3_nzero_global = 0;
   int64_t xray4_nzero_global = 0;
   int64_t temp_nzero_global = 0;
   int64_t rho_nzero_global = 0;
   int64_t phi_nzero_global = 0;
   int64_t vel_nzero_global = 0;


   double tsz_tot = 0;
   double ksz_tot = 0;
   double xray1_tot = 0;
   double xray2_tot = 0;
   double xray3_tot = 0;
   double xray4_tot = 0;
   double temp_tot = 0;
   double rho_tot = 0;
   double phi_tot = 0;
   double vel_tot = 0;

   double tsz_tot_global = 0;
   double ksz_tot_global = 0;
   double xray1_tot_global = 0;
   double xray2_tot_global = 0;
   double xray3_tot_global = 0;
   double xray4_tot_global = 0;
   double temp_tot_global = 0;
   double rho_tot_global = 0;
   double phi_tot_global = 0;
   double vel_tot_global = 0;

   double tsz_min = 1.e20;
   double ksz_min = 1.e20;
   double xray1_min = 1.e20;
   double xray2_min = 1.e20;
   double xray3_min = 1.e20;
   double xray4_min = 1.e20;
   double temp_min = 1.e20;
   double rho_min = 1.e20;
   double phi_min = 1.e20;
   double vel_min = 1.e20;
   

   double tsz_min_global = 1.e20;
   double ksz_min_global = 1.e20;
   double xray1_min_global = 1.e20;
   double xray2_min_global = 1.e20;
   double xray3_min_global = 1.e20;
   double xray4_min_global = 1.e20;
   double temp_min_global = 1.e20;
   double rho_min_global = 1.e20;
   double phi_min_global = 1.e20;
   double vel_min_global = 1.e20;


   double tsz_max = -1.e20;
   double ksz_max = -1.e20;
   double xray1_max = -1.e20;
   double xray2_max = -1.e20;
   double xray3_max = -1.e20;
   double xray4_max = -1.e20;
   double temp_max = -1.e20;
   double rho_max = -1.e20;
   double phi_max = -1.e20;
   double vel_max = -1.e20;

   double tsz_max_global = -1.e20;
   double ksz_max_global = -1.e20;
   double xray1_max_global = -1.e20;
   double xray2_max_global = -1.e20;
   double xray3_max_global = -1.e20;
   double xray4_max_global = -1.e20;
   double temp_max_global = -1.e20;
   double rho_max_global = -1.e20;
   double phi_max_global = -1.e20;
   double vel_max_global = -1.e20;


   for (int64_t j=0; j<ksz.size(); j++){
     tsz[j] = tsz[j]*TSZ_CONV;
     ksz[j] = ksz[j]*KSZ_CONV;
     xray1[j] = xray1[j]*XRAY_CONV;
     xray2[j] = xray2[j]*XRAY_CONV;
     xray3[j] = xray3[j]*XRAY_CONV;
     xray4[j] = xray4[j]*XRAY_CONV;

     tsz_tot += tsz[j];
     ksz_tot += ksz[j];
     xray1_tot += xray1[j];
     xray2_tot += xray2[j];
     xray3_tot += xray3[j];
     xray4_tot += xray4[j];
     temp_tot += temp[j];
     rho_tot += rho[j];
     phi_tot += phi[j];
     vel_tot += vel[j];

     if (tsz[j]!=0)
         tsz_nzero += 1;
     if (ksz[j]!=0)
         ksz_nzero += 1;
     if (xray1[j]!=0)
         xray1_nzero += 1;
     if (xray2[j]!=0)
         xray2_nzero += 1;
     if (xray3[j]!=0)
         xray3_nzero += 1;
     if (xray4[j]!=0)
         xray4_nzero += 1;
     if (temp[j]!=0)
         temp_nzero += 1;
     if (rho[j]!=0)
         rho_nzero += 1;
     if (phi[j]!=0)
         phi_nzero += 1;
     if (vel[j]!=0)
         vel_nzero += 1;
 
     tsz_min = (tsz[j]<tsz_min)?tsz[j]:tsz_min;
     ksz_min = (ksz[j]<ksz_min)?ksz[j]:ksz_min;
     xray1_min = (xray1[j]<xray1_min)?xray1[j]:xray1_min;
     xray2_min = (xray2[j]<xray2_min)?xray2[j]:xray2_min;
     xray3_min = (xray3[j]<xray3_min)?xray3[j]:xray3_min;
     xray4_min = (xray4[j]<xray4_min)?xray4[j]:xray4_min;
     temp_min = (ksz[j]<ksz_min)?ksz[j]:ksz_min;
     rho_min = (rho[j]<rho_min)?rho[j]:rho_min;
     phi_min = (phi[j]<phi_min)?phi[j]:phi_min;
     vel_min = (vel[j]<vel_min)?vel[j]:vel_min;

     tsz_max = (tsz[j]>tsz_max)?tsz[j]:tsz_max;
     ksz_max = (ksz[j]>ksz_max)?ksz[j]:ksz_max;
     xray1_max = (xray1[j]>xray1_max)?xray1[j]:xray1_max;
     xray2_max = (xray2[j]>xray2_max)?xray2[j]:xray2_max;
     xray3_max = (xray3[j]>xray3_max)?xray3[j]:xray3_max;
     xray4_max = (xray4[j]>xray4_max)?xray4[j]:xray4_max;
     temp_max = (ksz[j]>ksz_max)?ksz[j]:ksz_max;
     rho_max = (rho[j]>rho_max)?rho[j]:rho_max;
     phi_max = (phi[j]>phi_max)?phi[j]:phi_max;
     vel_max = (vel[j]>vel_max)?vel[j]:vel_max;

    } 
  // alter
  MPI_Reduce(&tsz_tot,&tsz_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&tsz_nzero,&tsz_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&tsz_min,&tsz_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&tsz_max,&tsz_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&ksz_tot,&ksz_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&ksz_nzero,&ksz_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&ksz_min,&ksz_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&ksz_max,&ksz_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&xray1_tot,&xray1_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray1_nzero,&xray1_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray1_min,&xray1_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray1_max,&xray1_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&xray2_tot,&xray2_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray2_nzero,&xray2_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray2_min,&xray2_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray2_max,&xray2_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&xray3_tot,&xray3_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray3_nzero,&xray3_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray3_min,&xray3_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray3_max,&xray3_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&xray4_tot,&xray4_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray4_nzero,&xray4_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray4_min,&xray4_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&xray4_max,&xray4_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&temp_tot,&temp_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&temp_nzero,&temp_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&temp_min,&temp_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&temp_max,&temp_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&rho_tot,&rho_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&rho_nzero,&rho_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&rho_min,&rho_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&rho_max,&rho_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&phi_tot,&phi_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&phi_nzero,&phi_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&phi_min,&phi_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&phi_max,&phi_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  MPI_Reduce(&vel_tot,&vel_tot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&vel_nzero,&vel_nzero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&vel_min,&vel_min_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&vel_max,&vel_max_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);


  if (commRank==0){
    // add outputs 
   cout << "map tsz sum = "<< tsz_tot_global << ", min = " << tsz_min_global << ", max = "<< tsz_max_global  << ", number of zeros = " << tsz_nzero_global <<  endl;
   cout << "map ksz sum = "<< ksz_tot_global << ", min = " << ksz_min_global << ", max = "<< ksz_max_global  << ", number of zeros = " << ksz_nzero_global <<  endl;
   cout << "map xray1 sum = "<< xray1_tot_global << ", min = " << xray1_min_global << ", max = "<< xray1_max_global  << ", number of zeros = " << xray1_nzero_global <<  endl;
   cout << "map xray2 sum = "<< xray2_tot_global << ", min = " << xray2_min_global << ", max = "<< xray2_max_global  << ", number of zeros = " << xray2_nzero_global <<  endl;
   cout << "map xray3 sum = "<< xray3_tot_global << ", min = " << xray3_min_global << ", max = "<< xray3_max_global  << ", number of zeros = " << xray3_nzero_global <<  endl;
   cout << "map xray4 sum = "<< xray4_tot_global << ", min = " << xray4_min_global << ", max = "<< xray4_max_global  << ", number of zeros = " << xray4_nzero_global <<  endl;
   cout << "map temp sum = "<< temp_tot_global << ", min = " << temp_min_global << ", max = "<< temp_max_global  << ", number of zeros = " << temp_nzero_global <<  endl;
   cout << "map rho sum = "<< rho_tot_global << ", min = " << rho_min_global << ", max = "<< rho_max_global  << ", number of zeros = " << rho_nzero_global <<  endl;
   cout << "map phi sum = "<< phi_tot_global << ", min = " << phi_min_global << ", max = "<< phi_max_global  << ", number of zeros = " << phi_nzero_global <<  endl;
   cout << "map vel sum = "<< vel_tot_global << ", min = " << vel_min_global << ", max = "<< vel_max_global  << ", number of zeros = " << vel_nzero_global <<  endl;


   /*cout << "particle tsz sum = "<< tsz_tot2*TSZ_CONV << endl;
   cout << "particle tsz sum (adiabatic with read mu ) = "<< tsz_tot3*TSZ_CONV << endl;
   cout << "particle tsz sum (adiabatic with mu0) = "<< tsz_tot4*TSZ_CONV << endl;
   cout << "xray sum = "<< xray1_tot << endl;*/

}

   return 0;
}



/* NOTE: do not use  - this needs updating to the new CLOUDY method
 *
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
*/

void clear_pixel_hydro(int64_t start_idx, int rank_diff, vector<double> &rho , vector<double> &phi, vector<double> &vel, vector<double> &ksz, vector<double> &tsz, vector<double> &xray1, vector<double> &xray2, vector<double> &xray3, vector<double> &xray4, vector<double> &temp){
    int64_t len_pixel = pow(4,rank_diff);
    for (int64_t idx=start_idx; idx<(len_pixel+start_idx); idx++){
        rho[idx] = 0.0;
	phi[idx] = 0.0;
	vel[idx] = 0.0;
	ksz[idx] = 0.0;
	tsz[idx] = 0.0;
	xray1[idx] = 0.0;
	xray2[idx] = 0.0;
        xray3[idx] = 0.0;
        xray4[idx] = 0.0;
        temp[idx] = 0.0;
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



void initialize_pixel_hydro(int pix_val,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, vector<double> &rho , vector<double> &phi, vector<double> &vel, vector<double> &ksz, vector<double> &tsz, vector<double> &xray1, vector<double> &xray2, vector<double> &xray3, vector<double> &xray4, vector<double> &temp, int64_t &count, vector<int64_t> &start_idx, vector<int64_t> &end_idx, vector<int64_t> &pixnum_start, vector<int64_t> &pixnum_end, int rank_diff,  unordered_map<int64_t, int64_t>* ring_to_idx){
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
			xray1.push_back(0.0);
			xray2.push_back(0.0);
			xray3.push_back(0.0);
			xray4.push_back(0.0);
			temp.push_back(0.0);
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
