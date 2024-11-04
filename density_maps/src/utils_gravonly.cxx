
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

#include <list>
#include <cmath>


//#include "chealpix.h"
#include "GenericIO.h"

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


#include "pointing.h"
#include "healpix_base.h"
#include "vec3.h"

#include "PLParticles.h"

#include "utils_gravonly.h"
#include "pix_funcs_gravonly.h"
#ifdef HACC_NVRAM
#include "HaccIOUtil.h"
#endif


using namespace gio;
using namespace std;


#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


static ptrdiff_t drand48elmt(ptrdiff_t i) {
  return ptrdiff_t(drand48()*i);
}



void read_particles(PLParticles* P, string file_name) {

  GenericIO GIO(MPI_COMM_WORLD,file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();

  (*P).Resize(num_elems + GIO.requestedExtraSpace());

  for (int i=0; i<N_FLOATS; ++i)
    GIO.addVariable((const string)float_names[i], *((*P).float_data[i]), true);
  for (int i=0; i<N_INTS; ++i)
    GIO.addVariable((const string)int_names[i], *((*P).int_data[i]), true);
  for (int i=0; i<N_INT64S; ++i)
    GIO.addVariable((const string)int64_names[i], *((*P).int64_data[i]), true);
  for (int i=0; i<N_DOUBLES; ++i)
    GIO.addVariable((const string)double_names[i], *((*P).double_data[i]), true);
  for (int i=0; i<N_MASKS; ++i)
    GIO.addVariable((const string)mask_names[i], *((*P).mask_data[i]), true);

  GIO.readData();
  (*P).Resize(num_elems);
}

void output_downsampled_particles(PLParticles* P, float downsampling_rate, string file_name, string input_file_name){

  GenericIO GIO(MPI_COMM_WORLD,input_file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  size_t subsize = round(downsampling_rate*num_elems);

  double PhysOrigin[3];
  double PhysScale[3];
  GIO.readPhysOrigin(PhysOrigin);
  GIO.readPhysScale(PhysScale);

#ifdef HACC_NVRAM
  string IOKey = "maps_downsample_hydro";
  HACCIOFlush(IOKey);//Flush any previous writes
  GenericIO NewGIO(MPI_COMM_WORLD, modifyPathNVRAM(file_name));
  NewGIO.setPartition(HACCNodeIOInfo::getIOData().file_num);
#else
  GenericIO NewGIO(MPI_COMM_WORLD,file_name);
#endif
  NewGIO.setNumElems(subsize);

  for (int d=0; d < 3; ++d){
          NewGIO.setPhysOrigin(PhysOrigin[d], d);
          NewGIO.setPhysScale(PhysScale[d], d);
  }

  vector<int64_t> indices(num_elems);
  for (int64_t i = 0; i<num_elems; i++){
  indices[i] = i;
  }
  vector<vector<float>> var_data_floats(N_FLOATS);
  vector<vector<int>> var_data_ints(N_INTS);
  vector<vector<int64_t>> var_data_int64s(N_INT64S);
  vector<vector<double>> var_data_doubles(N_DOUBLES);
  vector<vector<MASK_T>> var_data_masks(N_MASKS);

  for (int i = 0; i < N_FLOATS; ++i){
      var_data_floats[i].resize(subsize + NewGIO.requestedExtraSpace()/sizeof(float));
      for (int64_t j=0; j<subsize; j++){
          var_data_floats[i][j] =  (*P->float_data[i])[indices[j]];
      }
      NewGIO.addVariable((const string)float_names[i], var_data_floats[i], true); // no positional info
  }
  for (int i = 0; i < N_INTS; ++i){
      var_data_ints[i].resize(subsize + NewGIO.requestedExtraSpace()/sizeof(int));
      for (int64_t j=0; j<subsize; j++){
          var_data_ints[i][j] =  (*P->int_data[i])[indices[j]];
      }
      NewGIO.addVariable((const string)int_names[i], var_data_ints[i], true); // no positional info
  }
  for (int i = 0; i < N_INT64S; ++i){
      var_data_int64s[i].resize(subsize + NewGIO.requestedExtraSpace()/sizeof(int64_t));
      for (int64_t j=0; j<subsize; j++){
          var_data_int64s[i][j] =  (*P->int64_data[i])[indices[j]];
      }
      NewGIO.addVariable((const string)int64_names[i], var_data_int64s[i], true); // no positional info
  }
  for (int i = 0; i < N_DOUBLES; ++i){
      var_data_doubles[i].resize(subsize + NewGIO.requestedExtraSpace()/sizeof(double));
      for (int64_t j=0; j<subsize; j++){
          var_data_doubles[i][j] =  (*P->double_data[i])[indices[j]];
      }
      NewGIO.addVariable((const string)double_names[i], var_data_doubles[i], true); // no positional info
  }
  for (int i = 0; i < N_MASKS; ++i){
      var_data_masks[i].resize(subsize + NewGIO.requestedExtraSpace()/sizeof(MASK_T));
      for (int64_t j=0; j<subsize; j++){
          var_data_masks[i][j] =  (*P->mask_data[i])[indices[j]];
      }
      NewGIO.addVariable((const string)mask_names[i], var_data_masks[i], true); // no positional info
  }

  NewGIO.write();
#ifdef HACC_NVRAM
  HACCIOTransfer(IOKey, file_name);
#endif
  MPI_Barrier(MPI_COMM_WORLD);

}



template< class T >
static void distributeProperty(std::vector<T>& vA, T* sendbuf, 
                              const int64_t* vOrder,
                              size_t totalToSend, size_t totalToRecv, 
                              const int *sendcounts,
                              const int *sdispls, MPI_Datatype sendtype,
                              const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
                              MPI_Comm comm)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < totalToSend; ++i)
        sendbuf[i] = vA[ vOrder[i] ];

    vA.resize(totalToRecv);

    MPI_Alltoallv(&sendbuf[0], sendcounts, sdispls, sendtype,
                  &vA[0], recvcounts, rdispls, recvtype, comm);
}


// update send_counts and others to pass by reference to avoid memory copy
void redistribute_particles(PLParticles* P, vector<int> send_counts, int numranks,T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff){
    vector<int> recv_counts(numranks,0);
    vector<int> recv_offset(numranks,0);
    vector<int> send_offset(numranks,0);
    MPI_Alltoall( &send_counts[0], 1 , MPI_INT, &recv_counts[0], 1, MPI_INT, MPI_COMM_WORLD); // check this 

    for (int k=1;k<numranks;k++){
        send_offset[k] = send_offset[k-1] + send_counts[k-1];
        recv_offset[k] = recv_offset[k-1] + recv_counts[k-1];
    }

    // check this 
    size_t tot_send = send_offset.back()+send_counts.back();
    vector<int64_t> id_tot(tot_send);

    compute_ranks_index( P, map_lores, map_hires, numranks,  send_offset, id_tot,  rank_diff);

    // then for each variable   - loop over variables 
    vector<int64_t>scratch(tot_send);
    size_t tot_recv = recv_offset.back()+recv_counts.back();

	// check typings here - do more beyond the floats 
    for (int i=0; i<N_FLOATS; i++){
          string var_name = float_names[i];
	  distributeProperty((*P->float_data[i]), (float*)(&scratch[0]), &id_tot[0], tot_send, tot_recv, &send_counts[0], &send_offset[0], MPI_FLOAT, &recv_counts[0], &recv_offset[0], MPI_FLOAT, MPI_COMM_WORLD);
    }
    for (int i=0; i<N_DOUBLES; i++){
          string var_name = double_names[i];
          distributeProperty((*P->double_data[i]), (double*)(&scratch[0]), &id_tot[0], tot_send, tot_recv, &send_counts[0], &send_offset[0], MPI_DOUBLE, &recv_counts[0], &recv_offset[0], MPI_DOUBLE, MPI_COMM_WORLD);
    }
    for (int i=0; i<N_INTS; i++){
          string var_name = int_names[i];
          distributeProperty((*P->int_data[i]), (int*)(&scratch[0]), &id_tot[0], tot_send, tot_recv, &send_counts[0], &send_offset[0], MPI_INT, &recv_counts[0], &recv_offset[0], MPI_INT, MPI_COMM_WORLD);
    }
    for (int i=0; i<N_INT64S; i++){
          string var_name = int64_names[i];
          distributeProperty((*P->int64_data[i]), (int64_t*)(&scratch[0]), &id_tot[0], tot_send, tot_recv, &send_counts[0], &send_offset[0], MPI_INT64_T, &recv_counts[0], &recv_offset[0], MPI_INT64_T, MPI_COMM_WORLD);
    }
    for (int i=0; i<N_MASKS; i++){
          string var_name = mask_names[i];
          distributeProperty((*P->mask_data[i]), (uint16_t*)(&scratch[0]), &id_tot[0], tot_send, tot_recv, &send_counts[0], &send_offset[0], MPI_UINT16_T, &recv_counts[0], &recv_offset[0], MPI_UINT16_T, MPI_COMM_WORLD);
    }
    
    (*P).Resize(tot_recv);
    
    }




void read_and_redistribute(string file_name, int numranks, PLParticles* P,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff,bool output_downsampled, float downsampling_rate, string file_name_output){

    int commRank, commRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

    double t1, t2, t3, t4, t5;
	
    int status;
    vector<int> send_count(numranks,0);

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    read_particles(P, file_name);

    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();
    if (commRank==0){
      printf( "Read time is %f\n", t2 - t1 );
    fflush(stdout);
    }


    if (output_downsampled){
        assert(downsampling_rate<1.0);
        assert(downsampling_rate>0.0);
        output_downsampled_particles( P, downsampling_rate, file_name_output, file_name);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t3 = MPI_Wtime();
    if (commRank==0){
      printf( "Output downsampled particle time is %f\n", t3 - t2 );
    fflush(stdout);
    }


    status = compute_ranks_count( P,  map_lores, map_hires, numranks, send_count, rank_diff);
    MPI_Barrier(MPI_COMM_WORLD);
    t4 = MPI_Wtime();
    if (commRank==0){
      printf( "Time to compute ranks for communication is %f\n", t4 - t3 );
    fflush(stdout);
    }


    redistribute_particles(P, send_count,numranks,map_lores, map_hires, rank_diff);
    MPI_Barrier(MPI_COMM_WORLD);
    t5 = MPI_Wtime();
    if (commRank==0){
      printf( "Time for redistribution is %f\n", t5 - t4 );
    fflush(stdout);
    }



    //read_particles(P, file_name);
    //status = compute_ranks_count( P,  map_lores, map_hires, numranks, send_count, rank_diff); 
    //redistribute_particles(P, send_count,numranks,map_lores, map_hires, rank_diff);

    return;

}


int output_file(int rank, MPI_File &fh, MPI_Request &req , vector<float> &rho, vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pixnum_start){
    // output pixel values to MPI file with handle fh
    for (int i=0;i<start_idx.size();i++){
    int start_tmp = start_idx[i];
    int off_tmp = end_idx[i] - start_idx[i];
    MPI_Offset offset = pixnum_start[i]*sizeof(float);
    MPI_File_seek(fh,offset,MPI_SEEK_SET);
    MPI_File_iwrite(fh,&rho[start_tmp],off_tmp,MPI_FLOAT, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
    return 0;
}

int output_file_double(int rank, MPI_File &fh, MPI_Request &req , vector<double> &rho, vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pixnum_start){
    for (int i=0;i<start_idx.size();i++){
    int start_tmp = start_idx[i];
    int off_tmp = end_idx[i] - start_idx[i];
    MPI_Offset offset = pixnum_start[i]*sizeof(double);
    MPI_File_seek(fh,offset,MPI_SEEK_SET);
    MPI_File_iwrite(fh,&rho[start_tmp],off_tmp,MPI_DOUBLE, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
    return 0;
}



void write_files(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<double> rho, vector<double> phi, vector<double> vel, int64_t npix_hires){
  int commRank, commRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

  string IOKey = "healpix_maps";
  outfile = outfile + stepnumber + ".gio";
  outfile = modifyPath(outfile, IOKey, stepnumber);

  size_t size_n = 0;
  for(size_t i=0; i < start_idx.size();i++)size_n += end_idx[i] - start_idx[i];
  assert(size_n < std::numeric_limits<int>::max());
  assert(size_n == rho.size());

  size_t idx = 0;
  vector<int64_t>data_idx(size_n,0);
  vector<double>spare(size_n,0.0);
  for(size_t i=0; i < start_idx.size();i++){
    size_t del = end_idx[i] - start_idx[i];
    size_t start = start_idx[i];
    size_t id = pix_nums_start[i];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t j=0; j < del; j++){
      data_idx.at(idx+j) = id + j;
    }
    idx += del;
    assert(idx <= size_n);
  }
  assert(idx == size_n);

#ifdef HACC_NVRAM
  HACCIOFlush(IOKey);//Flush any previous writes
  GenericIO GIO(MPI_COMM_WORLD, modifyPathNVRAM(outfile));
  GIO.setPartition(HACCNodeIOInfo::getIOData().file_num);
#else
  GenericIO GIO(MPI_COMM_WORLD,outfile);
#endif
  size_t size_pad = size_n + std::max(GIO.requestedExtraSpace()/sizeof(double), GIO.requestedExtraSpace()/sizeof(int64_t));
  rho.resize(size_pad); phi.resize(size_pad); vel.resize(size_pad); 
  data_idx.resize(size_pad);
  GIO.setNumElems(size_n);
  GIO.addVariable("rho", rho, GenericIO::VarHasExtraSpace);
  GIO.addVariable("phi", phi, GenericIO::VarHasExtraSpace);
  GIO.addVariable("vel", vel, GenericIO::VarHasExtraSpace);
  GIO.addVariable("idx", data_idx, GenericIO::VarHasExtraSpace);
  GIO.write();
#ifdef HACC_NVRAM
  HACCIOTransfer(IOKey, outfile);
#endif
  rho.resize(size_n); phi.resize(size_n); vel.resize(size_n); 
  data_idx.resize(size_n);
}

/*void write_files(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<double> rho, vector<double> phi, vector<double> vel, int64_t npix_hires){ 

  int commRank, commRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

  string output_name = outfile + stepnumber + "_dens.bin";
  string output_name_phi = outfile + stepnumber + "_phi.bin";
  string output_name_vel = outfile + stepnumber + "_vel.bin";


  const char *name_out = output_name.c_str();
  const char *name_out_phi = output_name_phi.c_str();
  const char *name_out_vel = output_name_vel.c_str();

  MPI_File fh, fh_phi, fh_vel;
  MPI_Request req, req_phi, req_vel;

  MPI_File_open(MPI_COMM_WORLD, name_out, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  MPI_File_open(MPI_COMM_WORLD, name_out_phi, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_phi);
  MPI_File_open(MPI_COMM_WORLD, name_out_vel, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_vel);

  // PL NOTE: make sure this is in the right place
  MPI_Offset filesize_double = npix_hires*sizeof(double);
  MPI_File_set_size(fh, filesize_double);
  MPI_File_set_size(fh_phi, filesize_double);
  MPI_File_set_size(fh_vel, filesize_double); 


  int status;
  status = output_file_double(commRank, fh, req, rho, start_idx, end_idx,  pix_nums_start);
  status = output_file_double(commRank, fh_phi, req_phi, phi, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_vel, req_vel, vel, start_idx, end_idx, pix_nums_start);

  MPI_File_close(&fh);
  MPI_File_close(&fh_phi);
  MPI_File_close(&fh_vel);

  MPI_Barrier(MPI_COMM_WORLD);
  }*/
