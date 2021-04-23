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

#include "MurmurHashNeutral2.cpp" 


/// Assumes ids and  are consistent but ordering may have changed
/// Redistributes aongst ranks and sorts the values for a one-to-one check of the changes in 
/// several halo outputs


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

Halos_test H_1;
Halos_test H_2;

bool comp_by_fof_dest(const halo_properties_test &a, const halo_properties_test &b) {
  return a.rank < b.rank;
}

bool comp_by_fof_id(const halo_properties_test &a, const halo_properties_test &b) {
  return a.fof_halo_tag < b.fof_halo_tag;
}

int compute_mean_float(vector<float> *val1, vector<float> *val2, int num_halos , string var_name, string var_name2, float lim){
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  double diff=0;
  double diff_frac=0;
  int n_tot;
  double frac_max;
  double mean;
  double stddev=0;
  double stddevq=0;
  double meanq=0;
  double meanq_tot;

  for (int i=0; i<num_halos; i++){
      diff += (double)(val1->at(i)-val2->at(i));
      meanq += (double)(val1->at(i));
      if (val1->at(i)!=0){
        double frac = (double)(fabs(val1->at(i)-val2->at(i))/fabs(val1->at(i)));
         diff_frac = (diff_frac<frac)?frac:diff_frac;
      }
   }
      MPI_Barrier(Partition::getComm());
      MPI_Allreduce(&diff, &mean, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
      MPI_Allreduce(&diff_frac, &frac_max, 1, MPI_DOUBLE, MPI_MAX,  Partition::getComm());
      MPI_Allreduce(&meanq, &meanq_tot, 1, MPI_DOUBLE, MPI_SUM, Partition::getComm());
      MPI_Allreduce(&num_halos, &n_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());
   mean = mean/n_tot;
   meanq_tot = meanq_tot/n_tot;

   double diffq;
   double diff_tmp;
   for (int i=0; i< num_halos; i++){
      diff_tmp = (double)(val1->at(i)-val2->at(i))-mean;
      diffq = (double)(val1->at(i))- meanq_tot;
      stddev += diff_tmp*diff_tmp/(n_tot-1);
      stddevq += diffq*diffq/(n_tot-1);
   }
   double stddev_tot;
   double stddevq_tot;
   MPI_Allreduce(&stddev, &stddev_tot, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   MPI_Allreduce(&stddevq, &stddevq_tot, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   stddev_tot = sqrt(stddev_tot);
   stddevq_tot = sqrt(stddevq_tot);

   bool print_out = true;
   if (rank==0){

     if ((frac_max<lim)||((fabs(stddev_tot/stddevq_tot)<lim)&&(fabs(mean/meanq_tot)<lim))) // no values change by more than a percent
       print_out=false;
     if (print_out){
     cout << " " << var_name << endl;
     cout << " " << var_name2 << endl;
     cout << " ______________________________________" <<endl;
     cout << " mean difference = " << mean << endl;
     cout << " maximum fractional difference = " << frac_max<< endl;
     cout << " standard deviation of difference = " << stddev_tot << endl;
     cout << " mean of quantity = " << meanq_tot << endl;
     cout << " standard deviation of quantity = " << stddevq_tot << endl;
     cout << endl;
     return 1;
     }
   }
   return 0;

}

int  compute_mean_std_dist(Halos_test H_1 , Halos_test H_2, float lim ){
  // compute the mean and std of the differences and make a histogram to look more closely
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = H_1.num_halos;

  int err=0;
  for (int i =0; i<N_HALO_FLOATS; i++){
    string var_name = float_var_names_test[i];
    string var_name2 = float_var_names_test2[i];
    err += compute_mean_float(H_1.float_data[i],H_2.float_data[i],count,var_name,var_name2,lim);
  }
  return err;
}


inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}


void read_halos(Halos_test &H0, string file_name, int file_opt) {
 // Read halo files into a buffer
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();

  H0.Resize(num_elems + GIO.requestedExtraSpace()); 

  GIO.addVariable("fof_halo_tag",   *(H0.fof_halo_tag),true);
  GIO.addVariable("fof_halo_count", *(H0.fof_halo_count), true);
  if (H0.has_sod)
    GIO.addVariable("sod_halo_count", *(H0.sod_halo_count), true);

 if (file_opt==1){
  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names_test[i], *(H0.float_data[i]), true);
   }
 else{
  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names_test2[i], *(H0.float_data[i]), true);
 }

  GIO.readData();
  H0.Resize(num_elems);
}  

int main( int argc, char** argv ) {
  MPI_Init( &argc, &argv );
  Partition::initialize();
  GenericIO::setNaturalDefaultPartition();

  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  string fof_file     = string(argv[1]);
  string fof_file2    = string(argv[2]);
  stringstream thresh{ argv[3] };
  float lim{};
  if (!(thresh >> lim))
	  lim = 0.01;

  // Create halo buffers
  H_1.Allocate();
  H_1.has_sod = true;
  H_2.Allocate();
  H_2.has_sod = true;
  H_1.Set_MPIType();
  H_2.Set_MPIType();

  // Reading halos
  read_halos(H_1, fof_file, 1);
  read_halos(H_2, fof_file2, 2);



  // determine destination ranks
  vector<halo_properties_test> fof_halo_send;
  vector<int> fof_halo_send_cnt(n_ranks,0);
  vector<halo_properties_test> fof_halo_send2;
  vector<int> fof_halo_send_cnt2(n_ranks,0);

  for (int i=0; i<H_1.num_halos; ++i) {
    halo_properties_test tmp = H_1.GetProperties(i);
    tmp.rank = tag_to_rank(tmp.fof_halo_tag, n_ranks);
    fof_halo_send.push_back(tmp);
    ++fof_halo_send_cnt[tmp.rank];
  }
  for (int i=0; i<H_2.num_halos; ++i) {
    halo_properties_test tmp = H_2.GetProperties(i);
    tmp.rank = tag_to_rank(tmp.fof_halo_tag, n_ranks);
    fof_halo_send2.push_back(tmp);
    ++fof_halo_send_cnt2[tmp.rank];
  }

  // sort by destination rank
  sort(fof_halo_send.begin(),fof_halo_send.end(), comp_by_fof_dest);
  sort(fof_halo_send2.begin(),fof_halo_send2.end(), comp_by_fof_dest);
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
  std::sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_fof_id);
  std::sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_fof_id);



  // write into buffers
  H_1.Resize(0);
  H_2.Resize(0);
  for (int i=0; i<fof_halo_recv_total; ++i) {
    halo_properties_test tmp =  fof_halo_recv[i];
    H_1.PushBack(tmp);
  }
  for (int i=0; i<fof_halo_recv_total2; ++i){
       halo_properties_test tmp = fof_halo_recv2[i];
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
       if (H_1.fof_halo_tag->at(i)!=H_2.fof_halo_tag->at(i))
         err+=1;
     }
   }
   int numh1 = H_1.num_halos;
   int numh2 = H_2.num_halos;
   int dn_1 = 0;
   int dn_2 = 0;

   if (skip_err&&(err>0)){
       // if you want to skip over missing particles: note only do this if the catalogs should mostly match up.
       err = 0;
       int i=0;
       while (i < numh1) {
          if (H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i)){
             i += 1;
           }
          else {
           bool not_found = true;
           for (int j=0;j<32;j++){
               if (H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i+j)){
                  for (int k=0; k<j; k++)
                      H_2.Erase(i+k);
                not_found = false;
                i+=1; // iterate once found
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

  // compute the error characteristics for the catalogs
  err = compute_mean_std_dist(H_1 ,H_2, lim);

    if ((rank==0)&&(err==0)){
      cout << " Results " << endl;
      cout << " ______________________________ " << endl;
      cout << endl;      
      cout << " Comparison test passed! " << endl;
      cout << " All variables within threshold of "  << lim << endl;
      cout << " Total number of non-matching halos = "<< ndiff_tot+ndiff_tot2 << endl;
      cout << endl;
      cout << " ______________________________ " << endl;
  }
  if ((rank==0)&&(err>0)){
      cout << " Results " << endl;
      cout << " ______________________________ " << endl;
      cout << endl;
      cout << " Comparison exceeded threshold of " << lim << " for " << err << " variables" << endl;
      cout << " out of a total of " <<  N_HALO_FLOATS << " variables " << endl;
      cout << " See above outputs for details  "<< endl;
      cout << " Total number of non-matching halos = "<< ndiff_tot+ndiff_tot2 << endl;
      cout << endl;
      cout << " ______________________________ " << endl;
  }



  MPI_Barrier(Partition::getComm());
  H_1.Deallocate();
  H_2.Deallocate();


  Partition::finalize();
  MPI_Finalize();
  return 0;
}
