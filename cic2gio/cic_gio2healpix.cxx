#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>

#include "GenericIO.h"

#include "BasicDefinition.h"

#include <sstream>

//TODO add link to basic definitions
// cosmotools/algorithms/halofinder/BasicDefinition.h


#include <healpix_base.h>

#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "fitshandle.h"
//#include "share_utils.h"
#include "string_utils.h"
#include "healpix_tables.h"
//#include "c_utils.h"
#include "math_utils.h"
#include "vec3.h"
#include "stdio.h"
#include "fitsio.h"

#include <fstream>
#include <vector>
#include <algorithm>

#define MASK_SPECIES 2
#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t

// uncomment for adiabatic simulations
#define HYBRID_SG



using namespace std;
using namespace gio;

int main(int argc, char *argv[]) {

  int root_process = 0;
  MPI_Init(&argc, &argv);

  int commRank, commRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

  int num_proc = commRanks;
  double t1, t2, t3; 
  t1 = MPI_Wtime(); 

  if(argc < 8) {
    fprintf(stderr,"USAGE: %s <mpiioName> <outfile_ksz> <outfile_tsz> <sample_rate> <nsteps> <start_step> <nside> \n", argv[0]);
    exit(-1);
  }
  
  if(commRank == root_process) printf("commRanks and commrank are %i %i \n", commRanks, commRank);

  char *mpiioName_base = argv[1];
  char *outfile_ksz = argv[2];
  char *outfile_tsz = argv[3];
 
  int nsteps;
  float samplerate;
  int start_step;
  int nside_in; 

  samplerate = atof(argv[4]);
  nsteps = atoi(argv[5]);
  start_step = atoi(argv[6]);
  nside_in = atoi(argv[7]);
  assert(nside_in>1);
  assert(nside_in<32768);


  vector<POSVEL_T> xx, yy, zz;
  vector<POSVEL_T> vx, vy, vz;
  vector<POSVEL_T> a; 


#ifdef HYBRID_SG
  const double mu0 = MU_ION;
#else
  vector<POSVEL_T> mu;
#endif
  vector<POSVEL_T> mass; 
  vector<POSVEL_T> uu;

  assert(sizeof(ID_T) == 8);
  vector<ID_T> id;

  size_t Np = 0;
  unsigned Method = GenericIO::FileIOPOSIX;
  const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
  if (EnvStr && string(EnvStr) == "1")
    Method = GenericIO::FileIOMPI;
 
  Healpix_Ordering_Scheme ring = RING;
  int64 nside= (int64) nside_in;

  T_Healpix_Base<int64> map; // this certainly exists
  map = Healpix_Base2 (nside, ring, SET_NSIDE);
   
  int order = Healpix_Base::nside2order(nside);
  int64 npix = map.Npix();

  // area in steradians
  float pixsize = (4.*3.141529/npix);

  arr<double> map_output_total_ksz(npix);
  arr<double> map_output_total_tsz(npix);
  vector<double> map_output_ksz;
  vector<double> map_output_tsz;

  map_output_ksz.resize(npix);
  map_output_tsz.resize(npix);
  for (int64 i=0;i<npix;i++){
      map_output_ksz[i] = 0.0;
      map_output_tsz[i] = 0.0;
  }

  for (int ii=0;ii<nsteps;ii++) { 

    char step[10*sizeof(char)];
    sprintf(step,"%d",start_step+ii);

    char mpiioName[150];
    strcpy(mpiioName,mpiioName_base);
    strcat(mpiioName,step);


    { // scope GIO

      GenericIO GIO(MPI_COMM_WORLD, mpiioName, Method);
      GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
      MPI_Barrier(MPI_COMM_WORLD);

      Np = GIO.readNumElems();
      printf("rank = %d Np = %lu \n", commRank, Np);

      xx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      yy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      zz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      a.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      id.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      vx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      vy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      vz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      uu.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      mass.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
#ifndef HYBRID_SG
      mu.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
#endif

      //  READ IN VARIABLES

      GIO.addVariable("x", xx, true);
      GIO.addVariable("y", yy, true);
      GIO.addVariable("z", zz, true);
      GIO.addVariable("vx", vx, true);
      GIO.addVariable("vy", vy, true);
      GIO.addVariable("vz", vz, true);
      GIO.addVariable("mass", mass, true);
      GIO.addVariable("a", a, true);
      GIO.addVariable("uu", uu, true);
#ifndef HYBRID_SG
      GIO.addVariable("mu", mu, true);
#endif
      GIO.addVariable("id", id, true);

      GIO.readData();
    } // destroy GIO prior to calling MPI_Finalize

    printf("rank = %d Np_final = %lu \n", commRank, Np);
    xx.resize(Np);
    yy.resize(Np);
    zz.resize(Np);
    vx.resize(Np);
    vy.resize(Np);
    vz.resize(Np);
    uu.resize(Np);
    mass.resize(Np);
#ifndef HYBRID_SG
    mu.resize(Np);
#endif
    a.resize(Np);

    vector<float> total_weights;
    total_weights.resize(Np);

    double amin = 1.e20;
    double amax = -1.e20;
    double mumin = 1.e20;
    double mumax = -1.e20;
    double mmin = 1.e20;
    double mmax = -1.e20;
    double umin = 1.e20;
    double umax = -1.e20;
    double dcmin = 1.e20;
    double dcmax = -1.e20;

    for (int i=0; i<Np; i++) {

      //TODO: check if wind should be included  
      //if (!isGas(mask[i])) continue; // only use gas particles
      if (!isInterGas(mask[i])) continue;

      double xd = xx[i];
      double yd = yy[i];
      double zd = zz[i];

      double dist_comov2 = xd*xd + yd*yd + zd*zd;
      double vxd_los = vx[i] *(xd/sqrt(dist_comov2));
      double vyd_los = vy[i] *(yd/sqrt(dist_comov2));
      double vzd_los = vz[i] *(zd/sqrt(dist_comov2));  

      double v_los = vxd_los+vyd_los+vzd_los;
      //double vxd = vx[i];
      //double vyd = vy[i];
      //double vzd = vz[i];
     // double H = H0 * (om * a[i]*a[i]*a[i] + ol + o) ...
     //TODO: velocities might not have a factor of h! Check this! 
     // also v = d(ra)/dt = dr/dt a + r da/dt 
     //  v/a = dr/dt + r H 
     //  dr/dt = v/a - r H 

      vec3 vec_val = vec3(xd,yd,zd);

#ifdef HYBRID_SG
      double mui = mu0;
#else
      double mui = mu[i];
#endif
      double mi = mass[i];
      double ui = uu[i];
      double aa = a[i];

      amin = std::min(amin, aa) ; amax = std::max(amax, aa);
      mumin = std::min(mumin, mui) ; mumax = std::max(mumax, mui);
      mmin = std::min(mmin, mi) ; mmax = std::max(mmax, mi);
      umin = std::min(umin, ui) ; umax = std::max(umax, ui);
      dcmin = std::min(dcmin, dist_comov2) ; dcmax = std::max(dcmax, dist_comov2);
       
      //TODO make sure this is right 
      // int64 -  
      int64 pix_num = map.vec2pix(vec_val);
      map_output_ksz[pix_num] += mi*v_los/dist_comov2/aa; // one factor of a cancels from v_los and dist_comov2
      map_output_tsz[pix_num] += mi*mui*ui/dist_comov2;   // a^2 factors cancel out in ui and dist_comov2  
    }    

    printf("Finished accumulating particles for rank %d\n", commRank);

    double amin_g, amax_g;
    double mumin_g, mumax_g;
    double mmin_g, mmax_g;
    double umin_g, umax_g;
    double dcmin_g, dcmax_g;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&amin, &amin_g, 1, MPI_DOUBLE, MPI_MIN, root_process, MPI_COMM_WORLD); MPI_Reduce(&amax, &amax_g, 1, MPI_DOUBLE, MPI_MAX, root_process, MPI_COMM_WORLD);
    MPI_Reduce(&mumin, &mumin_g, 1, MPI_DOUBLE, MPI_MIN, root_process, MPI_COMM_WORLD); MPI_Reduce(&mumax, &mumax_g, 1, MPI_DOUBLE, MPI_MAX, root_process, MPI_COMM_WORLD);
    MPI_Reduce(&mmin, &mmin_g, 1, MPI_DOUBLE, MPI_MIN, root_process, MPI_COMM_WORLD); MPI_Reduce(&mmax, &mmax_g, 1, MPI_DOUBLE, MPI_MAX, root_process, MPI_COMM_WORLD);
    MPI_Reduce(&umin, &umin_g, 1, MPI_DOUBLE, MPI_MIN, root_process, MPI_COMM_WORLD); MPI_Reduce(&umax, &umax_g, 1, MPI_DOUBLE, MPI_MAX, root_process, MPI_COMM_WORLD);
    MPI_Reduce(&dcmin, &dcmin_g, 1, MPI_DOUBLE, MPI_MIN, root_process, MPI_COMM_WORLD); MPI_Reduce(&dcmax, &dcmax_g, 1, MPI_DOUBLE, MPI_MAX, root_process, MPI_COMM_WORLD);


    int istep = start_step+ii;
    if(commRank == root_process) std::cout << " step: " << istep << " amin: " << amin << " amax: " << amax << std::endl; 
    if(commRank == root_process) std::cout << " step: " << istep << " mumin: " << mumin << " mumax: " << mumax << std::endl; 
    if(commRank == root_process) std::cout << " step: " << istep << " mmin: " << mmin << " mmax: " << mmax << std::endl; 
    if(commRank == root_process) std::cout << " step: " << istep << " umin: " << umin << " umax: " << umax << std::endl; 
    if(commRank == root_process) std::cout << " step: " << istep << " dcmin: " << dcmin << " dcmax: " << dcmax << std::endl; 

  } // nsteps for loop

  // Physical constants (CGS units)
  const double SIGMAT = 6.652e-25;
  const double MP     = 1.6737236e-24;
  const double ME     = 9.1093836e-28;
  const double CLIGHT = 2.99792458e10;
  const double MUE    = 1.14;
  const double GAMMA  = 5.0/3.0;
  const double YP     = 0.24;
  const double NHE    = 0.0; // assume helium is neutral 
  const double CHIE   = (1.0-YP*(1.0-NHE/4.0))/(1.0-YP/2.0);
  const double MSUN_TO_G = 1.9885e33;
  const double KM_TO_CM  = 1.0e5;
  const double MPC_TO_CM = 3.0856776e24;

  // kSZ and tSZ conversion factors. After applying, both are dimensionless with one factor of h to be included later 
  const float KSZ_CONV = (float)(-SIGMAT*CHIE/MUE/MP/CLIGHT)*(MSUN_TO_G*KM_TO_CM/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize);
  const float TSZ_CONV = (float)((GAMMA-1.0)*SIGMAT*CHIE/MUE/ME/CLIGHT/CLIGHT)*(MSUN_TO_G*KM_TO_CM*KM_TO_CM/MPC_TO_CM/MPC_TO_CM)*(1.0/samplerate)*(1.0/pixsize);
  if(commRank == root_process) {
    std::cout << " CHIE: " << CHIE << " SAMPLERATE: " << samplerate << std::endl;
    std::cout << " KSZ_CONV: " << KSZ_CONV << " TSZ_CONV: " << TSZ_CONV << std::endl;
  }
  for (int64 j=0; j<npix; j++) {
    map_output_ksz[j] *= KSZ_CONV;
    map_output_tsz[j] *= TSZ_CONV;
  }

   // Accumulate over all ranks
  int ncell = 16;
  int npix2 = npix/ncell;

  for (int64 j=0; j<npix; j++) {
    map_output_total_ksz[j] = 0.0;
    map_output_total_tsz[j] = 0.0;
  }

  t3 = MPI_Wtime();
  if(commRank == root_process) printf(" ELAPSED TIME FOR FIRST STAGE: %f\n", t3 - t1 );
  if(commRank == root_process) std::cout << " STARTING REDUCE WITH npix: " << npix << " npix2: " << npix2 << std::endl; 

  for (int i=0; i<ncell; i++) {
    MPI_Reduce(&map_output_ksz[npix2*i], &map_output_total_ksz[npix2*i], npix2, MPI_DOUBLE, MPI_SUM, root_process, MPI_COMM_WORLD);
    MPI_Reduce(&map_output_tsz[npix2*i], &map_output_total_tsz[npix2*i], npix2, MPI_DOUBLE, MPI_SUM, root_process, MPI_COMM_WORLD);
    double sum_count_ksz = 0.0;
    double sum_count_tsz = 0.0;
    for (int k=0;k<npix2;k++) {
      sum_count_ksz += map_output_total_ksz[i*npix2+k];
      sum_count_tsz += map_output_total_tsz[i*npix2+k];
    }
    if(commRank == root_process) std::cout << " ncell: " << i << " sum_count_ksz: " << sum_count_ksz << " sum_count_tsz: " << sum_count_tsz << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  double sum_count_ksz = 0.0;
  double sum_count_tsz = 0.0;
  for (int64 j=0; j<npix; j++) {
    sum_count_ksz += (map_output_total_ksz[j]*map_output_total_ksz[j]); 
    sum_count_tsz += (map_output_total_tsz[j]*map_output_total_tsz[j]); 
  }
  sum_count_ksz = sqrt(sum_count_ksz)/npix;
  sum_count_tsz = sqrt(sum_count_tsz)/npix;
  if(commRank == root_process) std::cout << " mean squared y: " << sum_count_tsz << " b: " << sum_count_ksz << std::endl; 

  arr<float> map_output_total_ksz_float(npix);
  arr<float> map_output_total_tsz_float(npix);
  for (int64 j=0; j<npix; j++) {
    map_output_total_ksz_float[j] = (float)map_output_total_ksz[j]; 
    map_output_total_tsz_float[j] = (float)map_output_total_tsz[j];
  }

  if(commRank == root_process){
    Healpix_Map<float> map_new_ksz(map_output_total_ksz_float,ring);
    Healpix_Map<float> map_new_tsz(map_output_total_tsz_float,ring);
    if(commRank == root_process) printf("about to write ksz and tsz files\n");
    write_Healpix_map_to_fits(outfile_ksz,map_new_ksz,planckType<float>());
    write_Healpix_map_to_fits(outfile_tsz,map_new_tsz,planckType<float>());
  }

  t2 = MPI_Wtime();
  if(commRank == root_process) printf( "Elapsed time is %f\n", t2 - t1 );

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

