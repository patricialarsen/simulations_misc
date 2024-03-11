#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>


#include "GenericIO.h"
#include "BasicDefinition.h"

#include <sstream>

#include "chealpix.h"


#include <healpix_base.h>
#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "fitshandle.h"
#include "string_utils.h"
#include "healpix_tables.h"
#include "math_utils.h"
#include "vec3.h"
#include "stdio.h"
#include "fitsio.h"


#include <fstream>
#include <vector>
#include <algorithm>


#define MASK_SPECIES 2
#define POSVEL_T float
#define MASK_T uint16_t

// uncomment for adiabatic simulations
#define HYBRID_SG

using namespace std;
using namespace gio;


struct sz_props{
    int64_t pix_num;
    double tsz;
    double ksz;
    int rank;
};

struct sz_output{
    int64_t pix_num;
    float tsz;
    float ksz;
};

MPI_Datatype set_SZ_MPI_Type()
{
    MPI_Datatype SZ_MPI_Type;
    MPI_Datatype type[4] = { MPI_INT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };
    int blocklen[4] = {1,1,1,1};
    MPI_Aint disp[4] = {offsetof(sz_props,pix_num),
                        offsetof(sz_props,tsz),
                        offsetof(sz_props,ksz),
                        offsetof(sz_props,rank)};
    MPI_Type_struct(4,blocklen,disp,type,&SZ_MPI_Type);
    MPI_Type_commit(&SZ_MPI_Type);
    return SZ_MPI_Type;
}

/*
void Set_MPIType(MPI_Datatype SZ_MPI_Type){
    MPI_Datatype type[4] = { MPI_INT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };
    int blocklen[4] = {1,1,1,1};

    sz_props sz;

    MPI_Aint base;
    MPI_Aint disp[4];

    MPI_Get_address(&sz, &base);
    MPI_Get_address(&sz.pix_num,  &disp[0]);
    MPI_Get_address(&sz.tsz,   &disp[1]);
    MPI_Get_address(&sz.ksz,    &disp[2]);
    MPI_Get_address(&sz.rank,    &disp[3]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;

    MPI_Type_struct(4,blocklen,disp,type,&SZ_MPI_Type);
    MPI_Type_commit(&SZ_MPI_Type);
}

void Set_MPIType_OUT(MPI_Datatype SZ_MPI_Type){
    MPI_Datatype type[3] = { MPI_INT64_T, MPI_FLOAT, MPI_FLOAT };
    int blocklen[3] = {1,1,1};

    sz_output sz;

    MPI_Aint base;
    MPI_Aint disp[3];

    MPI_Get_address(&sz, &base);
    MPI_Get_address(&sz.pix_num,  &disp[0]);
    MPI_Get_address(&sz.tsz,   &disp[1]);
    MPI_Get_address(&sz.ksz,    &disp[2]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; 

    MPI_Type_struct(3,blocklen,disp,type,&SZ_MPI_Type);
    MPI_Type_commit(&SZ_MPI_Type);
}
*/
MPI_Datatype set_SZ_OUT_MPI_Type()
{
    MPI_Datatype SZ_OUT_MPI_Type;
    MPI_Datatype type[3] = { MPI_INT64_T, MPI_FLOAT, MPI_FLOAT};
    int blocklen[3] = {1,1,1};
    MPI_Aint disp[3] = {offsetof(sz_output,pix_num),
                        offsetof(sz_output,tsz),
                        offsetof(sz_output,ksz)};
    MPI_Type_struct(3,blocklen,disp,type,&SZ_OUT_MPI_Type);
    MPI_Type_commit(&SZ_OUT_MPI_Type);
    return SZ_OUT_MPI_Type;
}



bool comp_by_pixnum(const sz_output &a, const sz_output &b){
  return a.pix_num < b.pix_num;
}


bool comp_by_rank(const sz_props &a, const sz_props &b){
  return a.rank < b.rank;
}

int num_to_rank(int64_t pix_num, int commRanks){
   return pix_num%commRanks; 
}


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

  // assert user chooses sensible map resolution 
  assert(nside_in>1);
  assert(nside_in<32768);


  vector<POSVEL_T> xx, yy, zz;
  vector<POSVEL_T> vx, vy, vz;
  vector<POSVEL_T> a; 
  vector<MASK_T> mask;
#ifdef HYBRID_SG
  const double mu0 = MU_ION;
#else
  vector<POSVEL_T> mu;
#endif
  vector<POSVEL_T> mass; 
  vector<POSVEL_T> uu;


  size_t Np = 0;
  unsigned Method = GenericIO::FileIOPOSIX;
  const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
  if (EnvStr && string(EnvStr) == "1")
    Method = GenericIO::FileIOMPI;
 

  Healpix_Ordering_Scheme ring = RING;
  int64_t nside= (int64_t) nside_in;
  int64_t npix = nside2npix64(nside);
  float pixsize = (4.*3.141529/npix);

  // create MPI datatype
  //
  MPI_Datatype SZ_OUT_MPI_Type = set_SZ_OUT_MPI_Type();
  MPI_Datatype SZ_MPI_Type = set_SZ_MPI_Type();

  //MPI_Datatype SZ_MPI_Type;
  //Set_MPIType(SZ_MPI_Type);
  //MPI_Datatype SZ_OUT_MPI_Type;
  //Set_MPIType_OUT(SZ_OUT_MPI_Type);

  // create map for pixel list on rank 
  map<int64_t,int64_t> pixel_rank;
  map<int64_t,int64_t> pixel_rank_inv;

  int tmp_rank;
  int pix_count = 0; 
  for (int64_t i=0; i< npix; i++){
    tmp_rank = num_to_rank(i,commRanks);
    if (tmp_rank==commRank){
        pixel_rank[i] = pix_count;
        pixel_rank_inv[pix_count] = i; 
        pix_count++; 
    }
  }

    vector<double> map_output_ksz, map_output_tsz;
    map_output_ksz.resize(pix_count);
    map_output_tsz.resize(pix_count);
    for (int64_t i=0; i<pix_count; i++){
        map_output_ksz[i] = 0;
        map_output_tsz[i] = 0;
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
      vx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      vy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      vz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      uu.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      mass.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
      mask.resize(Np + GIO.requestedExtraSpace()/sizeof(MASK_T));
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
      GIO.addVariable("mask",mask,true);

      GIO.readData();
    } // destroy GIO prior to calling MPI_Finalize

    //printf("rank = %d Np_final = %lu \n", commRank, Np);
    xx.resize(Np);
    yy.resize(Np);
    zz.resize(Np);
    vx.resize(Np);
    vy.resize(Np);
    vz.resize(Np);
    uu.resize(Np);
    mass.resize(Np);
    mask.resize(Np);

#ifndef HYBRID_SG
    mu.resize(Np);
#endif
    a.resize(Np);


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


    vector<int> send_count;
    send_count.resize(commRanks,0);

    vector<sz_props> sz_send;
    vector<sz_props> sz_recv;
    sz_send.resize(Np);


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
      double vec_val[3]; 
      vec_val[0] = (double) xd; vec_val[1] = (double) yd; vec_val[2] = (double) zd;

#ifdef HYBRID_SG
      double mui = mu0;
#else
      double mui = mu[i];
#endif
      double mi = mass[i];
      double ui = uu[i];
      double aa = a[i];

      // check for min/max values as a sanity check 
      amin = std::min(amin, aa) ; amax = std::max(amax, aa);
      mumin = std::min(mumin, mui) ; mumax = std::max(mumax, mui);
      mmin = std::min(mmin, mi) ; mmax = std::max(mmax, mi);
      umin = std::min(umin, ui) ; umax = std::max(umax, ui);
      dcmin = std::min(dcmin, dist_comov2) ; dcmax = std::max(dcmax, dist_comov2);
       
      int64_t pix_num_tmp; 
      double tsz_pix_tmp;
      double ksz_pix_tmp;
      vec2pix_ring64(nside, vec_val, &pix_num_tmp); 

      sz_send[i].pix_num = pix_num_tmp;
      sz_send[i].tsz = mi*mui*ui/dist_comov2;
      sz_send[i].ksz = mi*v_los/dist_comov2/aa;
      sz_send[i].rank = num_to_rank(pix_num_tmp,commRanks);
      ++send_count[sz_send[i].rank];
    }    

    std::sort(sz_send.begin(),sz_send.end(),comp_by_rank);


    vector<int> recv_count;
    recv_count.resize(commRanks,0);
    MPI_Alltoall(&send_count[0],1,MPI_INT,&recv_count[0],1,MPI_INT,MPI_COMM_WORLD);

    // compute offsets
    vector<int> send_offset;
    send_offset.resize(commRanks,0);
    vector<int> recv_offset;
    recv_offset.resize(commRanks,0);
    send_offset[0] = recv_offset[0] = 0;
    for (int i=1; i<commRanks; ++i) {
      send_offset[i] = send_offset[i-1] + send_count[i-1];
      recv_offset[i] = recv_offset[i-1] + recv_count[i-1];
    }

    // compute totals
    int64_t recv_total = 0;
    for (int i=0; i<commRanks; ++i) {
      recv_total += recv_count[i];
    }

    sz_recv.resize(recv_total);


    MPI_Barrier(MPI_COMM_WORLD);
    if (commRank == 0) 
       cout << "About to send data" << endl;


    MPI_Alltoallv(&sz_send[0],&send_count[0],&send_offset[0],SZ_MPI_Type,\
                   &sz_recv[0],&recv_count[0],&recv_offset[0],SZ_MPI_Type,MPI_COMM_WORLD);


    sz_send.resize(0); // might not be optimal to keep resizing this

    // sum values 
    int64_t Np2 = sz_recv.size();
    for (int64_t i=0; i < Np2; i++){
      int64_t rank_pix = pixel_rank[sz_recv[i].pix_num];
      map_output_ksz[rank_pix] += sz_recv[i].ksz;
      map_output_tsz[rank_pix] += sz_recv[i].tsz;
    }

    // then multiply by scaling factor and finally combine into one map and output 

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
  //const double CLIGHT = 2.99792458e10;
  const double MUE    = 1.14;
  //const double GAMMA  = 5.0/3.0;
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
  for (int64_t j=0; j<pix_count; j++) {
    map_output_ksz[j] *= KSZ_CONV;
    map_output_tsz[j] *= TSZ_CONV;
  }


  // cast to float and get the pixel number 
  vector<sz_output> map_output;
  map_output.resize(pix_count); 
  for (int64_t j=0; j<pix_count; j++) {
    map_output[j].ksz = (float)map_output_ksz[j];
    map_output[j].tsz = (float)map_output_tsz[j];
    map_output[j].pix_num = pixel_rank_inv[j];
  }


  vector<sz_output> map_output2;
  map_output2.resize(pix_count);
  for (int64_t j=0; j<pix_count; j++) {
    map_output2[j].ksz = 0.0;
    map_output2[j].tsz = 0.0;
    map_output2[j].pix_num = 0;
  }

  // remove the doubles 
  map_output_tsz.resize(0);
  map_output_ksz.resize(0);

  // compute offsets - need to make sure the recv_offset fits into integer type - this will break for 16384 currently 
  vector<int> recv_offset_arr;
  vector<int> pix_count_arr;
  pix_count_arr.resize(commRanks);
  recv_offset_arr.resize(commRanks); 

  printf("pix count (1) for rank %d is %d \n ", commRank, pix_count);

  MPI_Allgather(&pix_count, 1, MPI_INT, &pix_count_arr[0],1,MPI_INT,MPI_COMM_WORLD);
  //MPI_Alltoall(&pix_count,1,MPI_INT,&pix_count_arr[0],1,MPI_INT,MPI_COMM_WORLD);

  printf("pix count for rank %d is %d \n ", commRank, pix_count_arr[commRank]);
  recv_offset_arr[0] = 0; 
  for (int i=1; i<commRanks; ++i) {
      recv_offset_arr[i] = recv_offset_arr[i-1] + pix_count_arr[i-1];
  }
  printf("recv offset for rank %d is %d \n ", commRank, recv_offset_arr[commRank]);


  // communicate map to rank 0 
  //
  /*
  if (commRank==0){
  MPI_Gatherv(MPI_IN_PLACE,pix_count_arr[commRank], SZ_OUT_MPI_Type,\
                   &map_output[0],&pix_count_arr[0],&recv_offset_arr[0],SZ_OUT_MPI_Type,0, MPI_COMM_WORLD);
}
   else {
      MPI_Gatherv(&map_output[0],pix_count_arr[commRank], SZ_OUT_MPI_Type,\
                   NULL,0,0,SZ_OUT_MPI_Type,0, MPI_COMM_WORLD);
}
*/

MPI_Gatherv(&map_output[0],pix_count_arr[commRank],SZ_OUT_MPI_Type,
            &map_output2[0],&pix_count_arr[0],&recv_offset_arr[0],SZ_OUT_MPI_Type,0,MPI_COMM_WORLD);

//  MPI_Gatherv(&map_output[0],pix_count_arr[commRank], SZ_OUT_MPI_Type,\
                   &map_output[0],&pix_count_arr[0],&recv_offset_arr[0],SZ_OUT_MPI_Type,0, MPI_COMM_WORLD);

  std::sort(map_output2.begin(),map_output2.end(),comp_by_pixnum);

 // now just have to write out the map 

  arr<float> map_output_total_ksz_float(npix);
  arr<float> map_output_total_tsz_float(npix);
  for (int64_t j=0; j<npix; j++) {
    map_output_total_ksz_float[j] = map_output2[j].ksz; 
    map_output_total_tsz_float[j] = map_output2[j].tsz;
  }
  if(commRank == root_process) printf("about to write ksz and tsz files\n");

  Healpix_Map<float> map_new_ksz(map_output_total_ksz_float,ring);
  Healpix_Map<float> map_new_tsz(map_output_total_tsz_float,ring);
  if(commRank == root_process) printf("about to write ksz and tsz files\n");
    write_Healpix_map_to_fits(outfile_ksz,map_new_ksz,planckType<float>());
    write_Healpix_map_to_fits(outfile_tsz,map_new_tsz,planckType<float>());



/*

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
  */



  t2 = MPI_Wtime();
  if(commRank == root_process) printf( "Elapsed time is %f\n", t2 - t1 );


  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Type_free(&SZ_MPI_Type);
  MPI_Type_free(&SZ_OUT_MPI_Type);

  MPI_Finalize();

  return 0;
}

