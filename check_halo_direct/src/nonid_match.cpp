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

/// Assumes halos match up and compares halos between mass bins but ids do not have to match
/// Redistributes amongst ranks and sorts the values for a one-to-one check of the changes in 
/// several halo outputs


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}


bool comp_by_fof_mass(const halo_properties_test &a, const halo_properties_test &b) {
  return a.float_data[0] < b.float_data[0];
}
bool comp_by_fof_x(const halo_properties_test &a, const halo_properties_test &b) {
  return a.float_data[2] < b.float_data[2];
}
bool comp_by_fof_y(const halo_properties_test &a, const halo_properties_test &b) {
  return a.float_data[3] < b.float_data[3];
}
bool comp_by_fof_z(const halo_properties_test &a, const halo_properties_test &b) {
  return a.float_data[4] < b.float_data[4];
}


bool check_comp_halo(const halo_properties_test &a, const halo_properties_test &b){
   float diff_x = a.float_data[2] - b.float_data[2];
   float diff_y  = a.float_data[3] - b.float_data[3];
   float diff_z = a.float_data[4] - b.float_data[4];
   float diff_m = (a.float_data[0] - b.float_data[0])/a.float_data[0];
   return ((fabs(diff_x)<0.1)&&(fabs(diff_y)<0.1)&&(fabs(diff_z)<0.1)&&(diff_m<0.1));
}



int vec_to_rank(float x, float y, float z, float box_size){
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();
  int nbins = pow(n_ranks,1./3);
  float bw = box_size/nbins;

  // correct boundaries 
  int xbin = (int)(x/bw); 
  int ybin = (int)(y/bw);
  int zbin = (int)(z/bw);
  if (xbin>2)
	  xbin = 2;
  if (xbin<0)
	  xbin=0;
  if (ybin>2)
	  ybin = 2;
  if (ybin<0)
	  ybin = 0;
  if (zbin>2)
	  zbin = 2;
  if (zbin<0)
	  zbin = 0;

  int dest_rank = xbin*nbins*nbins + ybin*nbins + zbin;
  return dest_rank;
}

int match_pos (string fof_file, string fof_file2, float lim, float box_size, float min_mass, float max_mass){
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  if (rank==0){
	  cout << " Comparing halo properties based on mass and position match "<< endl;
	  cout << " ___________________________________________________________ " << endl;
	  cout << endl;
  }

  Halos_test H_1;
  Halos_test H_2;

  // Create halo buffers
  H_1.Allocate();
  H_1.has_sod = true;
  H_2.Allocate();
  H_2.has_sod = true;
  H_1.Set_MPIType();
  H_2.Set_MPIType();


  // Reading halos
  read_halos(H_1, fof_file, 2);
  read_halos(H_2, fof_file2, 2);


  if (rank==0)
	  cout << endl;

  // determine destination ranks
  vector<halo_properties_test> fof_halo_send;
  vector<int> fof_halo_send_cnt(n_ranks,0);
  vector<halo_properties_test> fof_halo_send2;
  vector<int> fof_halo_send_cnt2(n_ranks,0);

  for (int64_t i=0; i<H_1.num_halos; ++i) {
    halo_properties_test tmp = H_1.GetProperties(i);
    tmp.rank = vec_to_rank(tmp.float_data[2],tmp.float_data[3],tmp.float_data[4],box_size);
    fof_halo_send.push_back(tmp);
    ++fof_halo_send_cnt[tmp.rank];
  }
  for (int64_t i=0; i<H_2.num_halos; ++i) {
    halo_properties_test tmp = H_2.GetProperties(i);
    tmp.rank = vec_to_rank(tmp.float_data[2],tmp.float_data[3],tmp.float_data[4],box_size);
    fof_halo_send2.push_back(tmp);
    ++fof_halo_send_cnt2[tmp.rank];
  }

  // sort by destination rank
  sort(fof_halo_send.begin(),fof_halo_send.end(), comp_by_halo_dest);
  sort(fof_halo_send2.begin(),fof_halo_send2.end(), comp_by_halo_dest);
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

  vector<halo_properties_test> fof_halo_recv, fof_halo_recv2;
  fof_halo_recv.resize(fof_halo_recv_total);
  fof_halo_recv2.resize(fof_halo_recv_total2);
  MPI_Barrier(Partition::getComm());
  //send data 
  MPI_Alltoallv(&fof_halo_send[0],&fof_halo_send_cnt[0],&fof_halo_send_off[0], H_1.halo_properties_MPI_Type,\
                   &fof_halo_recv[0],&fof_halo_recv_cnt[0],&fof_halo_recv_off[0], H_1.halo_properties_MPI_Type, Partition::getComm());

  MPI_Alltoallv(&fof_halo_send2[0],&fof_halo_send_cnt2[0],&fof_halo_send_off2[0], H_2.halo_properties_MPI_Type,\
                   &fof_halo_recv2[0],&fof_halo_recv_cnt2[0],&fof_halo_recv_off2[0], H_2.halo_properties_MPI_Type, Partition::getComm());



  // sort by fof halo tag
  //
  std::sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_fof_x);
  std::sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_fof_x);
  std::stable_sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_fof_y);
  std::stable_sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_fof_y);
  std::stable_sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_fof_z);
  std::stable_sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_fof_z);
  std::stable_sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_fof_mass);
  std::stable_sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_fof_mass);


  // write into buffers
  H_1.Resize(0);
  H_2.Resize(0);
  for (int64_t i=0; i<fof_halo_recv_total; ++i) {
    halo_properties_test tmp =  fof_halo_recv[i];
    if ((tmp.float_data[0]>min_mass)&&(tmp.float_data[0]<max_mass))
      H_1.PushBack(tmp);
  }
  for (int64_t i=0; i<fof_halo_recv_total2; ++i){
    halo_properties_test tmp = fof_halo_recv2[i];
    if ((tmp.float_data[0]>min_mass)&&(tmp.float_data[0]<max_mass))
      H_2.PushBack(tmp);
  }
  fof_halo_recv.resize(0);
  fof_halo_recv2.resize(0);


  int err = 1; // we do not assume that these match up 
  bool skip_err = true;

  if (H_1.num_halos != H_2.num_halos){
      err += 1;
    }

  int64_t numh1 = H_1.num_halos;
   int64_t numh2 = H_2.num_halos;
   int64_t numh1_i = H_1.num_halos;
   int dn_1 = 0;
   int dn_2 = 0;

   if (skip_err&&(err>0)){
       // if you want to skip over missing particles: note only do this if the catalogs should mostly match up.
       err = 0;
       int64_t i=0;
       while (i < numh1) {
	   halo_properties_test tmp = H_1.GetProperties(i);
	   halo_properties_test tmp2 = H_2.GetProperties(i);
          if (check_comp_halo(tmp,tmp2)){
             i += 1;
           }
          else {
           bool not_found = true;
           for (int j=0;j<10;j++){
	       if ((i+j)<H_2.num_halos){   
	       halo_properties_test tmp3 = H_2.GetProperties(i+j);
	       if (check_comp_halo(tmp,tmp3)){
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


   dn_1 = (int)(numh1_i- H_1.num_halos);
   dn_2 =  (int)(numh2 - H_2.num_halos);

  }
  

  int ndiff_tot = 0;
  int ndiff_tot2 = 0;
  int nh_tot = 0;
  MPI_Allreduce(&dn_1, &ndiff_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());
  MPI_Allreduce(&dn_2, &ndiff_tot2, 1, MPI_INT, MPI_SUM,  Partition::getComm());
  MPI_Allreduce(&numh1, &nh_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());


  // compute the error characteristics for the catalogs
  err = compute_mean_std_dist_halo(H_1 ,H_2, lim);

    if ((rank==0)&&(err==0)){
      cout << " Results " << endl;
      cout << " _______ " << endl;
      cout << endl;      
      cout << " Comparison test passed! " << endl;
      cout << " All variables within threshold of "  << lim << endl;
      cout << " Total number of non-matching halos = "<< ndiff_tot+ndiff_tot2 << endl;
      cout << " Total number of halos = "<<numh1 <<endl;
      cout << endl;
  }
  if ((rank==0)&&(err>0)){
      cout << " Results " << endl;
      cout << " _______ " << endl;
      cout << endl;
      cout << " Comparison exceeded threshold of " << lim << " for " << err << " variables" << endl;
      cout << " out of a total of " <<  N_HALO_FLOATS << " variables " << endl;
      cout << " See above outputs for details  "<< endl;
      cout << " Total number of non-matching halos = "<< ndiff_tot+ndiff_tot2 << endl;
      cout << " Total number of halos = "<<nh_tot <<endl;
      cout << endl;
  }



  MPI_Barrier(Partition::getComm());
  H_1.Deallocate();
  H_2.Deallocate();


  return 0;
}
