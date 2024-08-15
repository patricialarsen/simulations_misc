
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


#include "chealpix.h"
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

//#include "particle_def.h"
#include "PLParticles.h"
//#include "PLHalos.h"

#include "utils.h"
#include "pix_funcs.h"

using namespace gio;
using namespace std;


#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t

static ptrdiff_t drand48elmt(ptrdiff_t i) {
  return ptrdiff_t(drand48()*i);
}


void read_particles(PLParticles* P, string file_name) {

  GenericIO GIO(MPI_COMM_WORLD,file_name);
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

  GenericIO NewGIO(MPI_COMM_WORLD,file_name);
  NewGIO.setNumElems(subsize);

  for (int d=0; d < 3; ++d){
	  NewGIO.setPhysOrigin(PhysOrigin[d], d);
	  NewGIO.setPhysScale(PhysScale[d], d);
  }

  vector<int64_t> indices(num_elems);
  for (int64_t i = 0; i<num_elems; i++){
  indices[i] = i;
  }

  random_shuffle(indices.begin(), indices.end(), drand48elmt);

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




void read_and_redistribute(string file_name, int numranks, PLParticles* P,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff, bool output_downsampled, float downsampling_rate, string file_name_output){

    int status;
    vector<int> send_count(numranks,0);

    read_particles(P, file_name);
    
    if (output_downsampled){
        assert(downsampling_rate<1.0);
	assert(downsampling_rate>0.0);
	output_downsampled_particles( P, downsampling_rate, file_name_output, file_name);
    }

    status = compute_ranks_count( P,  map_lores, map_hires, numranks, send_count, rank_diff); 
    redistribute_particles(P, send_count,numranks,map_lores, map_hires, rank_diff);

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


void write_files_hydro(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<double> rho, vector<double> phi, vector<double> vel, vector<double> ksz, vector<double> tsz, vector<double> xray1, vector<double> xray2, vector<double> xray3, vector<double> xray4, vector<double> temp, int64_t npix_hires){

  int commRank, commRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

  string output_name = outfile + stepnumber + "_dens.bin";
  string output_name_phi = outfile + stepnumber + "_phi.bin";
  string output_name_vel = outfile + stepnumber + "_vel.bin";
  string output_name_ksz = outfile + stepnumber + "_ksz.bin";
  string output_name_tsz = outfile + stepnumber + "_tsz.bin";
  string output_name_x1 = outfile + stepnumber + "_xray_ROSAT.bin";
  string output_name_x2 = outfile + stepnumber + "_xray_ErositaLo.bin";
  string output_name_x3 = outfile + stepnumber + "_xray_ErositaHi.bin";
  string output_name_x4 = outfile + stepnumber + "_xray_Bolo.bin";
  string output_name_temp = outfile + stepnumber + "_temp.bin";


  const char *name_out = output_name.c_str();
  const char *name_out_phi = output_name_phi.c_str();
  const char *name_out_vel = output_name_vel.c_str();
  const char *name_out_ksz = output_name_ksz.c_str();
  const char *name_out_tsz = output_name_tsz.c_str();
  const char *name_out_x1 = output_name_x1.c_str();
  const char *name_out_x2 = output_name_x2.c_str();
  const char *name_out_x3 = output_name_x3.c_str();
  const char *name_out_x4 = output_name_x4.c_str();
  const char *name_out_temp = output_name_temp.c_str();


  MPI_File fh, fh_phi, fh_vel, fh_ksz, fh_tsz, fh_x1, fh_x2, fh_x3, fh_x4, fh_temp;
  MPI_Request req, req_phi, req_vel, req_ksz, req_tsz, req_x1, req_x2, req_x3, req_x4, req_temp;

  MPI_File_open(MPI_COMM_WORLD, name_out, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  MPI_File_open(MPI_COMM_WORLD, name_out_phi, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_phi);
  MPI_File_open(MPI_COMM_WORLD, name_out_vel, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_vel);
  MPI_File_open(MPI_COMM_WORLD, name_out_ksz, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_ksz);
  MPI_File_open(MPI_COMM_WORLD, name_out_tsz, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_tsz);
  MPI_File_open(MPI_COMM_WORLD, name_out_x1, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_x1);
  MPI_File_open(MPI_COMM_WORLD, name_out_x2, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_x2);
  MPI_File_open(MPI_COMM_WORLD, name_out_x3, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_x3);
  MPI_File_open(MPI_COMM_WORLD, name_out_x4, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_x4);
  MPI_File_open(MPI_COMM_WORLD, name_out_temp, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_temp);

  MPI_Offset filesize = npix_hires*sizeof(float);
  MPI_Offset filesize_double = npix_hires*sizeof(double);

  MPI_File_set_size(fh, filesize_double);
  MPI_File_set_size(fh_phi, filesize_double);
  MPI_File_set_size(fh_vel, filesize_double);
  MPI_File_set_size(fh_ksz, filesize_double);
  MPI_File_set_size(fh_tsz, filesize_double);
  MPI_File_set_size(fh_x1, filesize_double);
  MPI_File_set_size(fh_x2, filesize_double);
  MPI_File_set_size(fh_x3, filesize_double);
  MPI_File_set_size(fh_x4, filesize_double);
  MPI_File_set_size(fh_temp, filesize_double);


  int status;
  status = output_file_double(commRank, fh, req, rho, start_idx, end_idx,  pix_nums_start);
  status = output_file_double(commRank, fh_phi, req_phi, phi, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_vel, req_vel, vel, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_ksz, req_ksz, ksz, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_tsz, req_tsz, tsz, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_x1, req_x1, xray1, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_x2, req_x2, xray2, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_x3, req_x3, xray3, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_x4, req_x4, xray4, start_idx, end_idx, pix_nums_start);
  status = output_file_double(commRank, fh_temp, req_temp, temp, start_idx, end_idx, pix_nums_start);


  MPI_File_close(&fh);
  MPI_File_close(&fh_phi);
  MPI_File_close(&fh_vel);
  MPI_File_close(&fh_ksz);
  MPI_File_close(&fh_tsz);
  MPI_File_close(&fh_x1);
  MPI_File_close(&fh_x2);
  MPI_File_close(&fh_x3);
  MPI_File_close(&fh_x4);
  MPI_File_close(&fh_temp);

  MPI_Barrier(MPI_COMM_WORLD);
  }



void write_files(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<float> rho, vector<float> phi, vector<float> vel, int64_t npix_hires){ 

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
  MPI_Offset filesize = npix_hires*sizeof(float);
  MPI_File_set_size(fh, filesize);
  MPI_File_set_size(fh_phi, filesize);
  MPI_File_set_size(fh_vel, filesize); 


  int status;
  status = output_file(commRank, fh, req, rho, start_idx, end_idx,  pix_nums_start);
  status = output_file(commRank, fh_phi, req_phi, phi, start_idx, end_idx, pix_nums_start);
  status = output_file(commRank, fh_vel, req_vel, vel, start_idx, end_idx, pix_nums_start);

  MPI_File_close(&fh);
  MPI_File_close(&fh_phi);
  MPI_File_close(&fh_vel);

  MPI_Barrier(MPI_COMM_WORLD);
  }
