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

#include "Halos_test.h"

#include "MurmurHashNeutral2.h" 
#include "routines.h"

// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}


bool comp_by_part_dest(const particles_test &a, const particles_test &b) {
  return a.rank < b.rank;
}

bool comp_by_part_tag(const particles_test &a, const particles_test &b) {
  return a.fof_halo_tag < b.fof_halo_tag;
}
bool comp_by_part_id(const particles_test &a, const particles_test &b) {
  return a.id < b.id;
}



void read_particles(Particles_test &H0, string file_name) {
 // Read halo files into a buffer
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();

  H0.Resize(num_elems + GIO.requestedExtraSpace()); 

  GIO.addVariable("fof_halo_tag",   *(H0.fof_halo_tag),true);
  GIO.addVariable("id", *(H0.id), true);

  GIO.readData();
  H0.Resize(num_elems);
}  

int part_check(string fof_file, string fof_file2){

  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  if (rank==0){
      cout << endl;
      cout << " Checking constistency of halo tags and particle ids" << endl;
      cout << " ___________________________________________________" << endl;
      cout << endl; 
  }

  Particles_test H_1;
  Particles_test H_2;

  // Create halo buffers
  H_1.Allocate();
  H_2.Allocate();
  H_1.Set_MPIType();
  H_2.Set_MPIType();

  // Reading halos
  read_particles(H_1, fof_file);
  read_particles(H_2, fof_file2);

  if (rank==0)
	  cout << endl;

  // determine destination ranks
  vector<particles_test> fof_halo_send;
  vector<int> fof_halo_send_cnt(n_ranks,0);
  vector<particles_test> fof_halo_send2;
  vector<int> fof_halo_send_cnt2(n_ranks,0);

  for (int i=0; i<H_1.num_halos; ++i) {
    particles_test tmp = H_1.GetProperties(i);
    tmp.rank = tag_to_rank(tmp.fof_halo_tag, n_ranks);
    fof_halo_send.push_back(tmp);
    ++fof_halo_send_cnt[tmp.rank];
  }
  for (int i=0; i<H_2.num_halos; ++i) {
    particles_test tmp = H_2.GetProperties(i);
    tmp.rank = tag_to_rank(tmp.fof_halo_tag, n_ranks);
    fof_halo_send2.push_back(tmp);
    ++fof_halo_send_cnt2[tmp.rank];
  }

  // sort by destination rank
  sort(fof_halo_send.begin(),fof_halo_send.end(), comp_by_part_dest);
  sort(fof_halo_send2.begin(),fof_halo_send2.end(), comp_by_part_dest);
  MPI_Barrier(Partition::getComm());



  H_1.Resize(0);
  H_2.Resize(0);

  // create send and receive buffers and offsets
  vector<int> fof_halo_recv_cnt;
  vector<int> fof_halo_recv_cnt2;
  fof_halo_recv_cnt.resize(n_ranks,0);
  fof_halo_recv_cnt2.resize(n_ranks,0);

  MPI_Alltoall(&fof_halo_send_cnt[0],1,MPI_INT,&fof_halo_recv_cnt[0],1,MPI_INT,Partition::getComm());
  MPI_Alltoall(&fof_halo_send_cnt2[0],1,MPI_INT,&fof_halo_recv_cnt2[0],1,MPI_INT,Partition::getComm());

  vector<int> fof_halo_send_off;
  fof_halo_send_off.resize(n_ranks,0);
  vector<int> fof_halo_recv_off;
  fof_halo_recv_off.resize(n_ranks,0);
  vector<int> fof_halo_send_off2;
  fof_halo_send_off2.resize(n_ranks,0);
  vector<int> fof_halo_recv_off2;
  fof_halo_recv_off2.resize(n_ranks,0);


  fof_halo_send_off[0] = fof_halo_recv_off[0] = 0;
  fof_halo_send_off2[0] = fof_halo_recv_off2[0] = 0;

  for (int i=1; i<n_ranks; ++i) {
    fof_halo_send_off[i] = fof_halo_send_off[i-1] + fof_halo_send_cnt[i-1];
    fof_halo_recv_off[i] = fof_halo_recv_off[i-1] + fof_halo_recv_cnt[i-1];
    fof_halo_send_off2[i] = fof_halo_send_off2[i-1] + fof_halo_send_cnt2[i-1];
    fof_halo_recv_off2[i] = fof_halo_recv_off2[i-1] + fof_halo_recv_cnt2[i-1];
  }

  int fof_halo_recv_total = 0;
  int fof_halo_recv_total2 = 0;
  for (int i=0; i<n_ranks; ++i) {
    fof_halo_recv_total += fof_halo_recv_cnt[i];
    fof_halo_recv_total2 += fof_halo_recv_cnt2[i];
  }

  vector<particles_test> fof_halo_recv, fof_halo_recv2;
  fof_halo_recv.resize(fof_halo_recv_total);
  fof_halo_recv2.resize(fof_halo_recv_total2);
  MPI_Barrier(Partition::getComm());



  //send data 
  MPI_Alltoallv(&fof_halo_send[0],&fof_halo_send_cnt[0],&fof_halo_send_off[0], H_1.particles_test_MPI_Type,\
                   &fof_halo_recv[0],&fof_halo_recv_cnt[0],&fof_halo_recv_off[0], H_1.particles_test_MPI_Type, Partition::getComm());

  MPI_Alltoallv(&fof_halo_send2[0],&fof_halo_send_cnt2[0],&fof_halo_send_off2[0], H_2.particles_test_MPI_Type,\
                   &fof_halo_recv2[0],&fof_halo_recv_cnt2[0],&fof_halo_recv_off2[0], H_2.particles_test_MPI_Type, Partition::getComm());


  // sort by fof halo tag
  std::sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_part_tag);
  std::sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_part_tag);
  std::stable_sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_part_id);
  std::stable_sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_part_id);



  // write into buffers
  H_1.Resize(0);
  H_2.Resize(0);
  for (int i=0; i<fof_halo_recv_total; ++i) {
    particles_test tmp =  fof_halo_recv[i];
    H_1.PushBack(tmp);
  }
  for (int i=0; i<fof_halo_recv_total2; ++i){
       particles_test tmp = fof_halo_recv2[i];
    H_2.PushBack(tmp);
  }
  fof_halo_recv.resize(0);
  fof_halo_recv2.resize(0);


  int err = 0;
  bool skip_err = true;

  if (fof_halo_recv_total != fof_halo_recv_total2){
      err += 1;
    }

   else{
    for (int i=0;i<fof_halo_recv_total;i++){
       if ((H_1.fof_halo_tag->at(i)!=H_2.fof_halo_tag->at(i))||(H_1.id->at(i)!=H_2.id->at(i)))
         err+=1;
     }
   }

   int numh1 = H_1.num_halos;
   int numh2 = H_2.num_halos;
   int dn_1 = 0;
   int dn_2 = 0;
   if (skip_err&&(err>0)){
       err = 0;
       int i=0;
       while (i < numh1) {
          if ((H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i))&&(H_1.id->at(i)==H_2.id->at(i))){
             i += 1;
           }
          else {
           bool not_found = true;
           for (int j=0;j<32;j++){
               if((i+j)<H_2.num_halos){
               if ((H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i+j))&&(H_1.id->at(i)==H_2.id->at(i+j))){
                  for (int k=0; k<j; k++)
                      H_2.Erase(i+k);
                not_found = false;
                i+=1; // iterate once found
                }
	       }
              }
           if (not_found){
             H_1.Erase(i);
	     --numh1;
	   }
         }
     }

   dn_1 =  fof_halo_recv_total - H_1.num_halos;
   dn_2 =  fof_halo_recv_total2 - H_2.num_halos;

  }
  int ndiff_tot = 0;
  int ndiff_tot2 = 0;

  MPI_Allreduce(&dn_1, &ndiff_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());
  MPI_Allreduce(&dn_2, &ndiff_tot2, 1, MPI_INT, MPI_SUM,  Partition::getComm());

  if (rank==0){
  cout << " Results " << endl;
  cout << " _______ " << endl;
  cout << endl;      
  cout << " Total number of non-matching particles = "<< ndiff_tot+ndiff_tot2 << endl;
  cout << endl;
  }
  MPI_Barrier(Partition::getComm());
  H_1.Deallocate();
  H_2.Deallocate();


  return 0;
}
