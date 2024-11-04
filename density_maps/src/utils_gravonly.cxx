
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

int check_file(string filename){
    int commRank, commRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

    int valid = 0;
    if((commRank == 0) && (access(filename.c_str(), R_OK) == 0)){//using short-circuit
      valid = 1;
    }
    MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD); // step is valid
    return valid;
}


void read_particles(PLParticles* P, string file_name, string file_name_next) {

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


  if (check_file(file_name_next)){

  double t1 = MPI_Wtime();
  int rank, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  GenericIO GIO2(MPI_COMM_WORLD,file_name_next);
  GIO2.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems2 = GIO2.readNumElems();
  vector<int>rep2;
  vector<int64_t>id2;
  rep2.resize(num_elems2 + GIO2.requestedExtraSpace()/sizeof(int));
  id2.resize(num_elems2 + GIO2.requestedExtraSpace()/sizeof(int64_t));
  GIO2.addVariable("replication", rep2, true);
  GIO2.addVariable("id", id2, true);
  GIO2.readData();
  rep2.resize(num_elems2);
  id2.resize(num_elems2);

  /*add nicks stuff*/
  //Make a struct and MPI datatype to send out data
  struct Particle_dc {
    int64_t id;//id
    int rep;//replication
    int indx;//index in original array
    int rank;//rank we are sending to
    int home;//rank that it originated from
  };
  const int nitems = 5;
  int blocklengths[5] = {1, 1, 1, 1, 1};
  MPI_Aint offsets[5];
  MPI_Datatype types[5] = {MPI_INT64_T, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  struct Particle_dc particle;
  MPI_Aint base_address;
  MPI_Get_address(&particle, &base_address);
  MPI_Get_address(&particle.id, &offsets[0]);
  MPI_Get_address(&particle.rep, &offsets[1]);
  MPI_Get_address(&particle.indx, &offsets[2]);
  MPI_Get_address(&particle.rank, &offsets[3]);
  MPI_Get_address(&particle.home, &offsets[4]);
  for (int i = 0; i < nitems; i++) {
      offsets[i] -= base_address;
  }
  MPI_Datatype particle_type_dc;
  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &particle_type_dc);
  MPI_Type_commit(&particle_type_dc);

    //Map to rank based on id. Looks complicated im just jumbling id (using Thomas Wang's 64-bit integer hash function) so its balanced.
  auto id_to_rank = [](int64_t particle_id, int numRanks) -> int {
    uint64_t key = static_cast<uint64_t>(particle_id);
    key = (~key) + (key << 21);
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return static_cast<int>(key % numRanks);
  };

  vector<Particle_dc>P1_send;
  vector<Particle_dc>P1_recv;
  vector<int> P1_send_count(nRanks, 0);
  vector<int> P1_recv_count(nRanks, 0);
  vector<int> P1_send_offset(nRanks, 0);
  vector<int> P1_recv_offset(nRanks, 0);
  P1_send.resize(num_elems);
  assert(num_elems < numeric_limits<int>::max());
  assert(num_elems2 < numeric_limits<int>::max());
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(size_t i=0; i < num_elems; i++){
    int64_t id = (*P).int64_data[0]->at(i);// make sure index and rep are always first defined
    int rep = (*P).int_data[0]->at(i);
    int sendRank = id_to_rank(id,nRanks);
    P1_send[i] = {id, rep, int(i), sendRank, rank};
#ifdef _OPENMP
#pragma omp atomic
#endif
    P1_send_count[sendRank]++;
  }
  vector<Particle_dc>P2_send;
  vector<Particle_dc>P2_recv;
  vector<int> P2_send_count(nRanks, 0);
  vector<int> P2_recv_count(nRanks, 0);
  vector<int> P2_send_offset(nRanks, 0);
  vector<int> P2_recv_offset(nRanks, 0);
  P2_send.resize(num_elems2);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(size_t i=0; i < num_elems2; i++){
    int64_t id = id2.at(i);
    int rep = rep2.at(i);
    int sendRank = id_to_rank(id,nRanks);
    P2_send[i] = {id, rep, int(i), sendRank, rank};
#ifdef _OPENMP
#pragma omp atomic
#endif
    P2_send_count[sendRank]++;
  }
  //sort data by rank to send
  std::sort(P1_send.begin(), P1_send.end(), [](const Particle_dc& a, const Particle_dc& b) { return a.rank < b.rank;});
  std::sort(P2_send.begin(), P2_send.end(), [](const Particle_dc& a, const Particle_dc& b) { return a.rank < b.rank;});

  //Pass around sizes
  MPI_Alltoall( &P1_send_count[0], 1, MPI_INT, &P1_recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoall( &P2_send_count[0], 1, MPI_INT, &P2_recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);

  P1_send_offset[0]=0; P1_recv_offset[0]=0;
  P2_send_offset[0]=0; P2_recv_offset[0]=0;
  for (int i=1; i<nRanks; ++i) {
     P1_send_offset[i] = P1_send_offset[i-1] + P1_send_count[i-1];
     P1_recv_offset[i] = P1_recv_offset[i-1] + P1_recv_count[i-1];
     P2_send_offset[i] = P2_send_offset[i-1] + P2_send_count[i-1];
     P2_recv_offset[i] = P2_recv_offset[i-1] + P2_recv_count[i-1];
  }

  P1_recv.resize(P1_recv_offset.back() + P1_recv_count.back());
  P2_recv.resize(P2_recv_offset.back() + P2_recv_count.back());

  assert(P1_recv_offset.back() + P1_recv_count.back() < numeric_limits<int>::max());
  assert(P2_recv_offset.back() + P2_recv_count.back() < numeric_limits<int>::max());

  //Pass around data
  MPI_Alltoallv(&P1_send[0],&P1_send_count[0],&P1_send_offset[0],particle_type_dc,
                &P1_recv[0],&P1_recv_count[0],&P1_recv_offset[0],particle_type_dc, MPI_COMM_WORLD);
  MPI_Alltoallv(&P2_send[0],&P2_send_count[0],&P2_send_offset[0],particle_type_dc,
                &P2_recv[0],&P2_recv_count[0],&P2_recv_offset[0],particle_type_dc, MPI_COMM_WORLD);

  //Sort by id and replication
  std::sort(P1_recv.begin(), P1_recv.end(), [](const Particle_dc& a, const Particle_dc& b) {return (a.id < b.id) || (a.id == b.id && a.rep < b.rep);});
  std::sort(P2_recv.begin(), P2_recv.end(), [](const Particle_dc& a, const Particle_dc& b) {return (a.id < b.id) || (a.id == b.id && a.rep < b.rep);});

  P1_send.clear(); // Clear P1 buffer to store non-duplicate particles from P1 set
  std::fill(P1_send_count.begin(), P1_send_count.end(), 0);
  int p1_idx = 0;
  int p2_idx = 0;
  int p1_size = int(P1_recv.size());
  int p2_size = int(P2_recv.size());
  while (p1_idx < p1_size && p2_idx < p2_size) {
      const Particle_dc& p1 = P1_recv[p1_idx];
      const Particle_dc& p2 = P2_recv[p2_idx];
      if (p1.id == p2.id && p1.rep == p2.rep) { // Duplicate found
          p1_idx++;
          p2_idx++;
      } else if (p1.id < p2.id || (p1.id == p2.id && p1.rep < p2.rep)) { // Non-duplicate in P1
          P1_send.push_back(p1); // Save non-duplicate particle from set 1
          P1_send_count[p1.home]++;
          p1_idx++;
      } else { // Non-duplicate in P2
          p2_idx++;
      }
  }
  // Add remaining non-duplicate elements from P1_recv
  for (int i = p1_idx; i < p1_size; i++){
       const Particle_dc& p1 = P1_recv[i];
       P1_send.push_back(p1);
       P1_send_count[p1.home]++;
  }
  //Sort data by rank they originated from
  std::sort(P1_send.begin(), P1_send.end(), [](const Particle_dc& a, const Particle_dc& b) { return a.home < b.home;});

  //Send back sizes and data of first particle set
  MPI_Alltoall( &P1_send_count[0], 1, MPI_INT, &P1_recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);
  P1_send_offset[0]=0; P1_recv_offset[0]=0;
  for (int i=1; i<nRanks; ++i) {
     P1_send_offset[i] = P1_send_offset[i-1] + P1_send_count[i-1];
     P1_recv_offset[i] = P1_recv_offset[i-1] + P1_recv_count[i-1];
  }

  P1_recv.resize(P1_recv_offset.back() + P1_recv_count.back());

  MPI_Alltoallv(&P1_send[0],&P1_send_count[0],&P1_send_offset[0],particle_type_dc,
                &P1_recv[0],&P1_recv_count[0],&P1_recv_offset[0],particle_type_dc, MPI_COMM_WORLD);

  vector<int>indices(P1_recv.size());//Indices of final set
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(size_t i=0; i < P1_recv.size(); i++)indices[i] = P1_recv[i].indx;

  std::sort(indices.begin(), indices.end());//sort indices
  int64_t nDup = num_elems - indices.size(); int64_t nSize = num_elems;
  int64_t totalDup = 0, totalSize = 0;
  MPI_Allreduce(&nDup, &totalDup, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nSize, &totalSize, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  double t2 = MPI_Wtime();
  if(rank == 0)std::cout << "Total duplicates removed = " << totalDup << " total size = " << totalSize << " frac = " << double(totalDup)/double(totalSize) << " time = " << t2-t1 << " (s)" << std::endl;
  (*P).Transform(indices);//transform data

  MPI_Type_free(&particle_type_dc);//remove MPI datatype
  }
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




void read_and_redistribute(string file_name, string file_name_next, int numranks, PLParticles* P,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff,bool output_downsampled, float downsampling_rate, string file_name_output){

    int commRank, commRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

    double t1, t2, t3, t4, t5;
	
    int status;
    vector<int> send_count(numranks,0);

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    read_particles(P, file_name, file_name_next);

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

