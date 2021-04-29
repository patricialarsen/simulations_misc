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

#include "LC_test.h"

#include "MurmurHashNeutral2.cpp" 


/// Assumes ids and replication values are consistent but ordering may have changed
/// Redistributes amongst ranks and sorts the values for a one-to-one check of the changes in 
/// each of the main lightcone outputs 


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

LC_test L_1;
LC_test L_2;


bool comp_by_lc_dest(const lc_properties_test &a, const lc_properties_test &b){
  return a.rank < b.rank;
}

bool comp_by_id(const lc_properties_test &a, const lc_properties_test &b){
  return a.id < b.id;
}
bool comp_by_rep(const lc_properties_test &a, const lc_properties_test &b){
  return a.replication < b.replication;
}


int  compute_mean_std_dist(vector<float> *val1 , vector<float> *val2 , int count ,string var_name, float lim){
  // compute the mean and std of the differences and make a histogram to look more closely
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
  double meanq = 0;
  double meanq_tot;

  for (int i=0; i<count; i++){
      diff += (double)(val1->at(i)-val2->at(i));
      meanq += (double) (val1->at(i));
      if (val1->at(i)!=0){
        double frac = (double)(fabs(val1->at(i)-val2->at(i))/fabs(val1->at(i)));
         diff_frac = (diff_frac<frac)?frac:diff_frac;
      }
   }

   MPI_Allreduce(&diff, &mean, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   MPI_Allreduce(&diff_frac, &frac_max, 1, MPI_DOUBLE, MPI_MAX, Partition::getComm());
   MPI_Allreduce(&meanq, &meanq_tot, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   MPI_Allreduce(&count, &n_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());  
   mean = mean/n_tot;   
   meanq_tot = meanq_tot/n_tot;

   for (int i=0; i< count; i++){
      stddev += (double) ((val1->at(i)-val2->at(i)-(float)mean)*(val1->at(i)-val2->at(i)-(float)mean)/(n_tot-1));
      stddevq += ((double) (val1->at(i)) -meanq_tot)*((double)(val1->at(i)-meanq_tot))/( n_tot-1);
   }
   double stddev_tot, stddevq_tot;
   MPI_Reduce(&stddev, &stddev_tot, 1, MPI_DOUBLE, MPI_SUM, 0, Partition::getComm());
   MPI_Reduce(&stddevq, &stddevq_tot, 1, MPI_DOUBLE, MPI_SUM, 0, Partition::getComm());

   stddev_tot = sqrt(stddev_tot);
   stddevq_tot = sqrt(stddevq_tot);

   bool print_out = true;
   if (rank==0){
     if ((frac_max<lim)||((fabs(stddev_tot/stddevq_tot)<lim)&&(fabs(mean/meanq_tot)<lim))) // no values change by more than a percent
       print_out=false;
     if (print_out){
     cout << " " << var_name << endl;
     cout << " ______________________________________" <<endl;
     cout << " mean difference = "<< mean << endl;
     cout << " maximum fractional difference = "<< frac_max<< endl;
     cout << " standard deviation of difference = " << stddev_tot << endl;
     cout << endl;
     return 1;
     }
   }

  return 0;
}


int  compute_mean(LC_test L_1 , LC_test L_2, float lim ){
  // compute the mean and std of the differences and make a histogram to look more closely
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = L_1.num_parts;

  int err=0;
  for (int i =0; i<N_LC_FLOATS; i++){
    string var_name = float_var_names_test[i];
    err += compute_mean_std_dist(L_1.float_data[i],L_2.float_data[i],count,var_name,lim);
  }
  return err;
}



inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}


void read_lc_file(LC_test &L0, string file_name) {
  // reads in relevant lightcone information from a LC particle file
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  
  L0.Resize(num_elems + GIO.requestedExtraSpace());
 

  GIO.addVariable("id", *(L0.id), true);
  GIO.addVariable("replication", *(L0.replication), true);
  
  for (int i=0; i<N_LC_FLOATS; ++i)
	  GIO.addVariable((const string) float_var_names_test[i], *(L0.float_data[i]),true);

  GIO.readData();

  L0.Resize(num_elems);
}


int main( int argc, char** argv ) {
  MPI_Init( &argc, &argv );
  Partition::initialize();
  GenericIO::setNaturalDefaultPartition();

  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  string lc_file     = string(argv[1]);
  string lc_file2    = string(argv[2]);
  stringstream thresh{ argv[3] };
  float lim{};
  if (!(thresh >> lim))
	  lim = 0.01;


  L_1.Allocate();
  L_2.Allocate();
  L_1.Set_MPIType();
  L_2.Set_MPIType();



  int err = 0;
  bool skip_err = true;

  // LC HALOS
  read_lc_file(L_1, lc_file);                  // read the file
  read_lc_file(L_2, lc_file2);                  // read the file

  if (rank == 0)
    cout << "Done reading LC" << endl;

  vector<lc_properties_test> lc_halo_send;
  vector<int> lc_halo_send_cnt;
  lc_halo_send_cnt.resize(n_ranks,0);
  for (size_t i=0;i<L_1.num_parts;++i) {     // pack a vector of structures
	lc_properties_test tmp = L_1.GetProperties(i);
	tmp.rank = tag_to_rank(tmp.id,n_ranks);
	lc_halo_send.push_back(tmp);
	++lc_halo_send_cnt[tmp.rank];
  }                                       // prepare to send by contiguous destinations
  vector<lc_properties_test> lc_halo_send2;
  vector<int> lc_halo_send_cnt2;
  lc_halo_send_cnt2.resize(n_ranks,0);
  for (size_t i=0;i<L_2.num_parts;++i) {     // pack a vector of structures
        lc_properties_test tmp = L_2.GetProperties(i);
        tmp.rank = tag_to_rank(tmp.id,n_ranks);
        lc_halo_send2.push_back(tmp);
        ++lc_halo_send_cnt2[tmp.rank];
  }                                       // prepare to send by contiguous destinations


  if (rank == 0)
    cout << "Done packing LC buffers" << endl;

  std::sort(lc_halo_send.begin(),lc_halo_send.end(),comp_by_lc_dest);
  std::sort(lc_halo_send2.begin(),lc_halo_send2.end(),comp_by_lc_dest);

  // END LC HALOS

  if (rank == 0)
    cout << "Done sorting LC halos" << endl;

  L_1.Resize(0);
  L_2.Resize(0);



  // get the receive counts
  vector<int> lc_halo_recv_cnt;
  lc_halo_recv_cnt.resize(n_ranks,0);
  MPI_Alltoall(&lc_halo_send_cnt[0],1,MPI_INT,&lc_halo_recv_cnt[0],1,MPI_INT,Partition::getComm());

  vector<int> lc_halo_recv_cnt2;
  lc_halo_recv_cnt2.resize(n_ranks,0);
  MPI_Alltoall(&lc_halo_send_cnt2[0],1,MPI_INT,&lc_halo_recv_cnt2[0],1,MPI_INT,Partition::getComm());
  // each rank now knows how many items it will receive from every other rank

  // calculate the offsets
  vector<int> lc_halo_send_off;
  lc_halo_send_off.resize(n_ranks,0);
  vector<int> lc_halo_send_off2;
  lc_halo_send_off2.resize(n_ranks,0);


  vector<int> lc_halo_recv_off;
  lc_halo_recv_off.resize(n_ranks,0);
  vector<int> lc_halo_recv_off2;
  lc_halo_recv_off2.resize(n_ranks,0);

  lc_halo_send_off[0] = lc_halo_recv_off[0] = 0;
  lc_halo_send_off2[0] = lc_halo_recv_off2[0] = 0;

  for (int i=1; i<n_ranks; ++i) {
    lc_halo_send_off[i] = lc_halo_send_off[i-1] + lc_halo_send_cnt[i-1];
    lc_halo_recv_off[i] = lc_halo_recv_off[i-1] + lc_halo_recv_cnt[i-1];
    lc_halo_send_off2[i] = lc_halo_send_off2[i-1] + lc_halo_send_cnt2[i-1];
    lc_halo_recv_off2[i] = lc_halo_recv_off2[i-1] + lc_halo_recv_cnt2[i-1];
  }

  // compute the receive totals to allocate buffers
  int lc_halo_recv_total = 0;
  int lc_halo_recv_total2 = 0;

  for (int i=0; i<n_ranks; ++i) {
    lc_halo_recv_total += lc_halo_recv_cnt[i];
    lc_halo_recv_total2 += lc_halo_recv_cnt2[i];

  }

  vector<lc_properties_test> lc_halo_recv;
  lc_halo_recv.resize(lc_halo_recv_total);
  vector<lc_properties_test> lc_halo_recv2;
  lc_halo_recv2.resize(lc_halo_recv_total2);

  MPI_Barrier(Partition::getComm());

  if (rank == 0)
    cout << "About to send data" << endl;


  // send the actual data
  MPI_Alltoallv(&lc_halo_send[0],&lc_halo_send_cnt[0],&lc_halo_send_off[0],L_1.LC_properties_MPI_Type,\
                   &lc_halo_recv[0],&lc_halo_recv_cnt[0],&lc_halo_recv_off[0],L_1.LC_properties_MPI_Type,Partition::getComm());

  MPI_Alltoallv(&lc_halo_send2[0],&lc_halo_send_cnt2[0],&lc_halo_send_off2[0],L_2.LC_properties_MPI_Type,\
                   &lc_halo_recv2[0],&lc_halo_recv_cnt2[0],&lc_halo_recv_off2[0],L_2.LC_properties_MPI_Type,Partition::getComm());

  if (rank == 0)
    cout << "About to sort" << endl;

   std::sort(lc_halo_recv.begin(),lc_halo_recv.end(),comp_by_id);
   std::stable_sort(lc_halo_recv.begin(),lc_halo_recv.end(),comp_by_rep);
   std::sort(lc_halo_recv2.begin(),lc_halo_recv2.end(),comp_by_id);
   std::stable_sort(lc_halo_recv2.begin(),lc_halo_recv2.end(),comp_by_rep);

   // at this point we have saved all the data 
  if (rank == 0)
    cout << "Sorted" << endl;

  L_1.Resize(0);
  L_2.Resize(0);

  for (int i=0; i<lc_halo_recv_total; ++i){
     lc_properties_test tmp = lc_halo_recv[i];
     L_1.PushBack(tmp);
  }
    for (int i=0; i<lc_halo_recv_total2; ++i){
     lc_properties_test tmp = lc_halo_recv2[i];
     L_2.PushBack(tmp);
  }
    lc_halo_recv.resize(0);
    lc_halo_recv2.resize(0);


  if (rank == 0)
    cout << "Assigned to buffers" << endl;


  if (lc_halo_recv_total != lc_halo_recv_total2){
    err += 1;
      cout << "The number of elements is different in the two files for rank " << rank << endl;
      cout << lc_halo_recv_total << " , and " << lc_halo_recv_total2 << " for rank "<< rank<< endl;
    }
   else{

    for (int i=0;i<lc_halo_recv_total;i++){
       if( (L_1.id->at(i)!=L_2.id->at(i) )|| (L_1.replication->at(i)!=L_2.replication->at(i) ) )
         err+=1;
     }
   }

    int numl1 = L_1.num_parts;
    int numl2 = L_2.num_parts;
    int dn_1 = 0;
    int dn_2 = 0;


    if (skip_err&&(err>0)){
       err = 0;
       cout << "Skipping over missing particles for rank "<< rank << endl;
       int i = 0;
       while (i < numl1){
          if ((L_1.id->at(i)==L_2.id->at(i))&&(L_1.replication->at(i)==L_2.replication->at(i))){
             i += 1;
           }
          else {
           bool not_found = true;
           for (int j=0;j<32;j++){
               if ((L_1.id->at(i)==L_2.id->at(i+j))&&(L_1.replication->at(i)==L_2.replication->at(i+j))){
                  for (int k=0; k<j; k++)
                      L_2.Erase(i+k);
                not_found = false;
                i+=1;
                }
              }
           if (not_found){
             L_1.Erase(i); 
	     -- numl1;
	   }
         }
     }
     dn_1 = lc_halo_recv_total - L_1.num_parts;
     dn_2 = lc_halo_recv_total2 - L_2.num_parts;

    }

    int ndiff_tot = 0;
    int ndiff_tot2 = 0;
    MPI_Allreduce(&dn_1, &ndiff_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());
    MPI_Allreduce(&dn_2, &ndiff_tot2, 1, MPI_INT, MPI_SUM,  Partition::getComm());



    err = compute_mean(L_1,L_2,lim);

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
      cout << " out of a total of " <<  N_LC_FLOATS << " variables " << endl;
      cout << " See above outputs for details  "<< endl;
      cout << " Total number of non-matching particles = "<< ndiff_tot+ndiff_tot2 << endl;
      cout << endl;
      cout << " ______________________________ " << endl;
  }


  
 
 L_1.Deallocate();
 L_2.Deallocate();

  MPI_Barrier(Partition::getComm());


  Partition::finalize();
  MPI_Finalize();
  return 0;
}
