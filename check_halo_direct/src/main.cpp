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

typedef struct fof_halo {
  float mass;
  int64_t id;
  unsigned int destination_rank;
} fof_halo;


typedef struct halo_struct {
  float posvel_a_m[7];
  float phi;
  int rr;
  int64_t id;
  unsigned int destination_rank;
} lc_halo;

bool comp_by_lc_dest(const lc_halo &a, const lc_halo &b) {
  return a.destination_rank < b.destination_rank;
}

bool comp_by_id(const lc_halo &a, const lc_halo &b) {
  return a.id < b.id;
}

bool comp_by_rep(const lc_halo &a, const lc_halo &b) {
  return a.rr < b.rr;
}

bool comp_by_fof_dest(const halo_properties_test &a, const halo_properties_test &b) {
  return a.rank < b.rank;
}

bool comp_by_fof_id(const halo_properties_test &a, const halo_properties_test &b) {
  return a.fof_halo_tag < b.fof_halo_tag;
}

void compute_mean_float(vector<float> *val1, vector<float> *val2, int num_halos , string var_name){
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  double diff=0;
  double diff_frac=0;
  int n_tot;
  double frac_max;
  double mean;
  double stddev;

  for (int i=0; i<num_halos; i++){
      diff += (double)(val1->at(i)-val2->at(i));
      if (val1->at(i)!=0){
        double frac = (double)(fabs(val1->at(i)-val2->at(i))/fabs(val1->at(i)));
         diff_frac = (diff_frac<frac)?frac:diff_frac;
      }
   }

      MPI_Reduce(&diff, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, Partition::getComm());
      MPI_Reduce(&diff_frac, &frac_max, 1, MPI_DOUBLE, MPI_MAX, 0, Partition::getComm());
      MPI_Reduce(&num_halos, &n_tot, 1, MPI_INT, MPI_SUM, 0, Partition::getComm());
   mean = mean/n_tot;

   for (int i=0; i< num_halos; i++){
      stddev += (double) ((val1->at(i)-val2->at(i)-(float)mean)*(val1->at(i)-val2->at(i)-(float)mean)/(n_tot-1));
   }
   double stddev_tot;
   MPI_Reduce(&stddev, &stddev_tot, 1, MPI_DOUBLE, MPI_SUM, 0, Partition::getComm());
   stddev_tot = sqrt(stddev_tot);

   if (rank==0){
     cout << var_name << endl;
     cout << "______________________________________" <<endl;
     cout << " mean difference = "<< mean << endl;
     cout << " maximum fractional difference = "<< frac_max<< endl;
     cout << " standard deviation of difference = " << stddev_tot << endl;
     cout << endl;
   }

}

void  compute_mean_std_dist(Halos_test H_1 , Halos_test H_2 ){
  // compute the mean and std of the differences and make a histogram to look more closely
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = H_1.num_halos;

//    this->float_data[i]->at(idx)
  for (int i =0; i<N_HALO_FLOATS; i++){
    string var_name = float_var_names_test[i];
    compute_mean_float(H_1.float_data[i],H_2.float_data[i],count,var_name);
  }

  return;
}


inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}


void read_halos(Halos_test &H0, string file_name, int file_opt) {
 // may need a pointer to H0
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();

  H0.Resize(num_elems + GIO.requestedExtraSpace()); // think this works, may need a pointer

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




  H_1.Allocate();
  H_1.has_sod = true;
  H_2.Allocate();
  H_2.has_sod = true;

  // LC HALOS
  read_halos(H_1, fof_file, 1);
  read_halos(H_2, fof_file2, 2);

  if (rank == 0)
    cout << "Done reading halos" << endl;



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

  sort(fof_halo_send.begin(),fof_halo_send.end(), comp_by_fof_dest);
  sort(fof_halo_send2.begin(),fof_halo_send2.end(), comp_by_fof_dest);


  if (rank == 0)
    cout << "Sorted halos " << endl;

  H_1.Resize(0);
  H_2.Resize(0);

  MPI_Barrier(Partition::getComm());
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


  if (rank == 0)
    cout << "About to send data" << endl;

 // MPI_Datatype halo_properties_MPI_Type;
 // get_MPI_Type(halo_properties_MPI_Type);

MPI_Datatype halo_properties_MPI_Type;
{
    MPI_Datatype type[5] = { MPI_INT64_T, MPI_INT, MPI_INT64_T, MPI_INT, MPI_FLOAT };
    int blocklen[5] = {1,1,1,1,N_HALO_FLOATS};
    halo_properties_test hp;

    MPI_Aint base;
    MPI_Aint disp[5];

    MPI_Get_address(&hp, &base);
    MPI_Get_address(&hp.fof_halo_tag,     &disp[0]);
    MPI_Get_address(&hp.fof_halo_count,   &disp[1]);
    MPI_Get_address(&hp.sod_halo_count,   &disp[2]);
    MPI_Get_address(&hp.rank,             &disp[3]);
    MPI_Get_address(&hp.float_data,       &disp[4]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;
    disp[4]-=base;


    MPI_Type_struct(5,blocklen,disp,type,&halo_properties_MPI_Type);
    MPI_Type_commit(&halo_properties_MPI_Type);
}


  MPI_Alltoallv(&fof_halo_send[0],&fof_halo_send_cnt[0],&fof_halo_send_off[0], halo_properties_MPI_Type,\
                   &fof_halo_recv[0],&fof_halo_recv_cnt[0],&fof_halo_recv_off[0], halo_properties_MPI_Type, Partition::getComm());

  MPI_Alltoallv(&fof_halo_send2[0],&fof_halo_send_cnt2[0],&fof_halo_send_off2[0], halo_properties_MPI_Type,\
                   &fof_halo_recv2[0],&fof_halo_recv_cnt2[0],&fof_halo_recv_off2[0], halo_properties_MPI_Type, Partition::getComm());

  if (rank == 0)
    cout << "About to sort" << endl;

  std::sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_fof_id);
  std::sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_fof_id);


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

  int err = 0;
  bool skip_err = true;


  if (fof_halo_recv_total != fof_halo_recv_total2){
      err += 1;
      cout << "The number of elements is different in the two files for rank " << rank << endl;
      cout << fof_halo_recv_total << " , and " << fof_halo_recv_total2 << " for rank "<< rank<< endl;
    }
   else{
    for (int i=0;i<fof_halo_recv_total;i++){
       if (H_1.fof_halo_tag->at(i)!=H_2.fof_halo_tag->at(i))
         err+=1;
     }
   }
   if (skip_err&&(err>0)){
   cout << "Skipping over missing particles for rank "<<rank << endl;
  }



   int numh1 = H_1.num_halos;
   int numh2 = H_2.num_halos;

   int min_n = min(numh1,numh2);


    if (skip_err&&(err>0)){
       err = 0;
       cout << "Skipping over missing particles for rank "<< rank << endl;

       int i=0;
       while (i < min_n) {
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
                i+=1;
                }
              }
           if (not_found)
             H_1.Erase(i);
         }
     }
    }

    cout<< "number of elements before was " <<  fof_halo_recv_total2 << " , and "<<fof_halo_recv_total <<endl;
    cout<< "number of elements now is " <<  H_1.num_halos << " , and "<< H_2.num_halos <<endl;


    fof_halo_recv_total = H_1.num_halos;
    fof_halo_recv_total2 = H_2.num_halos;

    compute_mean_std_dist(H_1 ,H_2);


  MPI_Barrier(Partition::getComm());
  H_1.Deallocate();
  H_2.Deallocate();


  Partition::finalize();
  MPI_Finalize();
  return 0;
}
