#include <cstdlib>
#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdexcept>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <algorithm>    // std::sort
#include <sstream>
#include <omp.h>
#include <sys/stat.h>

// Generic IO
#include "GenericIO.h"
#include "Partition.h"

#include "HaloDistribute_mt.h"
#include "Halos.h"

#include "MurmurHashNeutral2.cpp" 

// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

//NEW
Halos H0;
cosmotk::HaloDistribute HD(&H0);

typedef struct lc_halo {
  float posvel_a_m[7];
  int rr[1];
  int64_t id;
  unsigned int destination_rank;
} lc_halo;

typedef struct fof_halo {
  float mass;
  int64_t id;
  unsigned int destination_rank;
} fof_halo;

bool comp_by_lc_dest(const lc_halo &a, const lc_halo &b) {
  return a.destination_rank < b.destination_rank;
}

//bool comp_by_fof_dest(const fof_halo &a, const fof_halo &b) {
bool comp_by_fof_dest(const halo_properties_t &a, const halo_properties_t &b) {
  return a.rank < b.rank;
}


//bool comp_by_fof_id(const fof_halo &a, const fof_halo &b) {
bool comp_by_fof_id(const halo_properties_t &a, const halo_properties_t &b) {
  return a.fof_halo_tag < b.fof_halo_tag;
}

struct IO_Buffers { // this is the order of lc_halo struct above

  // LC halo data
  vector<float> x;
  vector<float> y;
  vector<float> z;
  vector<float> vx;
  vector<float> vy;
  vector<float> vz;
  vector<float> a;
  vector<int> replication;
  vector<int64_t> id;
  double box_size[3];
  double origin[3];

  // 
};

IO_Buffers IOB;

void clear_IO_buffers() {
  IOB.x.clear();
  IOB.y.clear();
  IOB.z.clear();
  IOB.vx.clear();
  IOB.vy.clear();
  IOB.vz.clear();
  IOB.a.clear();
  IOB.replication.clear();
  IOB.id.clear();
}

inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}

void read_lc_file(string file_name) {
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  GIO.readPhysScale(IOB.box_size);
  GIO.readPhysOrigin(IOB.origin);
  IOB.x.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.y.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.z.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.vx.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.vy.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.vz.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.a.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.replication.resize(num_elems + GIO.requestedExtraSpace()/sizeof(int));
  IOB.id.resize(num_elems + GIO.requestedExtraSpace()/sizeof(int64_t));
  GIO.addVariable("x", IOB.x, true);
  GIO.addVariable("y", IOB.y, true);
  GIO.addVariable("z", IOB.z, true);
  GIO.addVariable("vx", IOB.vx, true);
  GIO.addVariable("vy", IOB.vy, true);
  GIO.addVariable("vz", IOB.vz, true);
  GIO.addVariable("a", IOB.a, true);
  GIO.addVariable("replication", IOB.replication, true);
  GIO.addVariable("id", IOB.id, true);
  GIO.readData();
  IOB.x.resize(num_elems);
  IOB.y.resize(num_elems);
  IOB.z.resize(num_elems);
  IOB.vx.resize(num_elems);
  IOB.vy.resize(num_elems);
  IOB.vz.resize(num_elems);
  IOB.a.resize(num_elems);
  IOB.replication.resize(num_elems);
  IOB.id.resize(num_elems);
}


void write_lc_file(string file_name) {
  GenericIO GIO(Partition::getComm(), file_name);
  GIO.setNumElems(IOB.x.size());
  GIO.setPhysOrigin(IOB.origin[0]);
  GIO.setPhysScale(IOB.box_size[0]);
  GIO.addVariable("x",IOB.x);
  GIO.addVariable("y",IOB.y);
  GIO.addVariable("z",IOB.z);
  GIO.addVariable("vx",IOB.vx);
  GIO.addVariable("vy",IOB.vy);
  GIO.addVariable("vz",IOB.vz);
  GIO.addVariable("a",IOB.a);
  GIO.addVariable("replication",IOB.replication);


  GIO.addVariable("fof_halo_tag", *(H0.fof_halo_tag));
  //GIO.addVariable("fof_halo_count", *(H0.fof_halo_count));
  
  //if (H0.has_sod)
    //GIO.addVariable("sod_halo_count", *(H0.sod_halo_count));

  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names[i], *(H0.float_data[i]));

  GIO.write();
}


halo_properties_t fof_lookup(int64_t tag, vector<halo_properties_t> &fof_halo_recv) {
  halo_properties_t tmp;
  tmp.fof_halo_tag = tag;
  vector<halo_properties_t>::iterator item = lower_bound(fof_halo_recv.begin(),fof_halo_recv.end(), tmp, comp_by_fof_id);
  if (item!=fof_halo_recv.end() && !comp_by_fof_id(tmp, *item)){
    return *item;
    }
  else{
    cout << "warning - halo not found " << endl;
    return tmp; 
}
}

int main( int argc, char** argv ) {
  MPI_Init( &argc, &argv );
  Partition::initialize();
  GenericIO::setNaturalDefaultPartition();

  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  string lc_file     = string(argv[1]);
  string fof_file    = string(argv[2]);
  string lc_out_file = string(argv[3]);





#ifdef _OPENMP



   if (rank==0){
    int nthreads, tid;



#pragma omp parallel private(nthreads, tid)
 {

 /* Obtain thread number */
 tid = omp_get_thread_num();
 printf("Hello World from thread = %d\n", tid);

 /* Only master thread does this */
 if (tid == 0) 
 {
  nthreads = omp_get_num_threads();
 printf("Number of threads = %d\n", nthreads);
 }

 } /* All threads join master thread and disband */


     printf("OpenMP is defined \n");
   }
#endif



  //NEW
  //HD.set_parameters(redistribute, Program.pad_factor, Simulation.rl, Program.halo_loss_threshold);
  HD.set_parameters(true, 0.0, 0.0, 0.0);
  HD.initialize();

  H0.Allocate();
  H0.has_sod = true;


  // LC HALOS
  read_lc_file(lc_file);                  // read the file

  if (rank == 0)
    cout << "Done reading LC" << endl;

  vector<lc_halo> lc_halo_send;
  vector<int> lc_halo_send_cnt;
  lc_halo_send_cnt.resize(n_ranks,0);
  for (size_t i=0;i<IOB.id.size();++i) {     // pack a vector of structures
    int64_t tag = IOB.id[i];
    if (tag<0)
      tag = (-1*tag) & 0x0000ffffffffffff;

    unsigned int recv_rank = tag_to_rank(tag, n_ranks);

    lc_halo h = { {IOB.x[i],IOB.y[i],IOB.z[i],IOB.vx[i],IOB.vy[i],IOB.vz[i],IOB.a[i]}, \
                  {IOB.replication[i]},\
                  tag, recv_rank};
    lc_halo_send.push_back(h);
    ++lc_halo_send_cnt[recv_rank];        // count the send items
  }                                       // prepare to send by contiguous destinations

  if (rank == 0)
    cout << "Done packing LC halos" << endl;

  std::sort(lc_halo_send.begin(),lc_halo_send.end(),comp_by_lc_dest);
  // END LC HALOS

  if (rank == 0)
    cout << "Done sorting LC halos" << endl;

  clear_IO_buffers();
  MPI_Barrier(Partition::getComm());

  HD.read_halos(fof_file);
  if (rank == 0)
    cout << "Read fof halos" << endl;

  vector<halo_properties_t> fof_halo_send;
  vector<int> fof_halo_send_cnt(n_ranks,0);
  
  // pack a buffer of halos with tag_to_rank for destinations and redistribute
  for (int i=0; i<H0.num_halos; ++i) {
    halo_properties_t tmp = H0.GetProperties(i);

    if (tmp.fof_halo_tag<0)
      tmp.fof_halo_tag = (-1*tmp.fof_halo_tag) & 0x0000ffffffffffff;

    tmp.rank = tag_to_rank(tmp.fof_halo_tag, n_ranks);


    fof_halo_send.push_back(tmp);
    ++fof_halo_send_cnt[tmp.rank];
  }

  if (rank == 0)
    cout << "Packed buffer" << endl;


  sort(fof_halo_send.begin(),fof_halo_send.end(), comp_by_fof_dest);
  

  MPI_Barrier(Partition::getComm());

  // create the MPI types
  MPI_Datatype lc_halo_type;
  {
    MPI_Datatype type[5] = { MPI_FLOAT, MPI_INT, MPI_INT64_T, MPI_UNSIGNED, MPI_UB };
    int blocklen[5] = {7,1,1,1,1};
    MPI_Aint disp[5] = {  offsetof(lc_halo,posvel_a_m),
                          offsetof(lc_halo,rr),
                          offsetof(lc_halo,id),
                          offsetof(lc_halo,destination_rank),
                          sizeof(lc_halo) };
    MPI_Type_struct(5,blocklen,disp,type,&lc_halo_type);
    MPI_Type_commit(&lc_halo_type);
  }


  // get the receive counts
  vector<int> lc_halo_recv_cnt;
  lc_halo_recv_cnt.resize(n_ranks,0);
  MPI_Alltoall(&lc_halo_send_cnt[0],1,MPI_INT,&lc_halo_recv_cnt[0],1,MPI_INT,Partition::getComm());

  vector<int> fof_halo_recv_cnt;
  fof_halo_recv_cnt.resize(n_ranks,0);
  MPI_Alltoall(&fof_halo_send_cnt[0],1,MPI_INT,&fof_halo_recv_cnt[0],1,MPI_INT,Partition::getComm());

  // each rank now knows how many items it will receive from every other rank

  // calculate the offsets
  vector<int> lc_halo_send_off;
  lc_halo_send_off.resize(n_ranks,0);
  vector<int> fof_halo_send_off;
  fof_halo_send_off.resize(n_ranks,0);
  vector<int> lc_halo_recv_off;
  lc_halo_recv_off.resize(n_ranks,0);
  vector<int> fof_halo_recv_off;
  fof_halo_recv_off.resize(n_ranks,0);

  lc_halo_send_off[0] = lc_halo_recv_off[0] = 0;
  fof_halo_send_off[0] = fof_halo_recv_off[0] = 0;
  for (int i=1; i<n_ranks; ++i) {
    lc_halo_send_off[i] = lc_halo_send_off[i-1] + lc_halo_send_cnt[i-1];
    lc_halo_recv_off[i] = lc_halo_recv_off[i-1] + lc_halo_recv_cnt[i-1];
    fof_halo_send_off[i] = fof_halo_send_off[i-1] + fof_halo_send_cnt[i-1];
    fof_halo_recv_off[i] = fof_halo_recv_off[i-1] + fof_halo_recv_cnt[i-1];
  }

  // compute the receive totals to allocate buffers
  int lc_halo_recv_total = 0;
  int fof_halo_recv_total = 0;
  for (int i=0; i<n_ranks; ++i) {
    lc_halo_recv_total += lc_halo_recv_cnt[i];
    fof_halo_recv_total += fof_halo_recv_cnt[i];
  }

  vector<lc_halo> lc_halo_recv;
  lc_halo_recv.resize(lc_halo_recv_total);
  //vector<fof_halo> fof_halo_recv;
  vector<halo_properties_t> fof_halo_recv;
  fof_halo_recv.resize(fof_halo_recv_total);
  MPI_Barrier(Partition::getComm());

  if (rank == 0)
    cout << "About to send data" << endl;


  // send the actual data
  MPI_Alltoallv(&lc_halo_send[0],&lc_halo_send_cnt[0],&lc_halo_send_off[0],lc_halo_type,\
                   &lc_halo_recv[0],&lc_halo_recv_cnt[0],&lc_halo_recv_off[0],lc_halo_type,Partition::getComm());

  MPI_Alltoallv(&fof_halo_send[0],&fof_halo_send_cnt[0],&fof_halo_send_off[0], HD.halo_properties_MPI_Type,\
                   &fof_halo_recv[0],&fof_halo_recv_cnt[0],&fof_halo_recv_off[0], HD.halo_properties_MPI_Type, Partition::getComm());
  if (rank == 0)
    cout << "About to sort" << endl;

  std::sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_fof_id);


  clear_IO_buffers();
  H0.Resize(0);
  for (int i=0; i<lc_halo_recv_total; ++i) {



    halo_properties_t tmp = fof_lookup(lc_halo_recv[i].id, fof_halo_recv);
    H0.PushBack(tmp);
    IOB.x.push_back(lc_halo_recv[i].posvel_a_m[0]);
    IOB.y.push_back(lc_halo_recv[i].posvel_a_m[1]);
    IOB.z.push_back(lc_halo_recv[i].posvel_a_m[2]);
    IOB.vx.push_back(lc_halo_recv[i].posvel_a_m[3]);
    IOB.vy.push_back(lc_halo_recv[i].posvel_a_m[4]);
    IOB.vz.push_back(lc_halo_recv[i].posvel_a_m[5]);
    IOB.a.push_back(lc_halo_recv[i].posvel_a_m[6]);
    IOB.replication.push_back(lc_halo_recv[i].rr[0]);
  }

  if (rank == 0)
    cout << "Filled IO buffer" << endl;


  MPI_Barrier(Partition::getComm());
  write_lc_file(lc_out_file);
  if (rank == 0)
    cout << "Written LC file" << endl;

  // NEW
  H0.Deallocate();
  HD.finalize();

  Partition::finalize();
  MPI_Finalize();
  return 0;
}
