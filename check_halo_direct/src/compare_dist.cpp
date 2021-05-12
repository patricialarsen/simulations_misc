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


int compute_mean_float_dist(vector<float> *val1, vector<float> *val2, int num_halos , int num_halos2, string var_name, string var_name2, float lim){
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  double stddev=0;
  double stddev2=0;
  double diff=0;
  double diff2=0;
  double max_val=0;
  double max_val2=0;
  int n_tot, n_tot2;
  double max2,max1;
  double mean, mean2;
  double stddev_tot, stddev_tot2;
  int err = 0;

  // mean computation
  for (int i=0; i<num_halos; i++){
     double tmp = (double)(val1->at(i));
     diff += tmp;
     max_val =(max_val<fabs(tmp))?fabs(tmp):max_val;
   }
    for (int i=0; i<num_halos2; i++){
     double tmp = (double)(val2->at(i));
     diff2 += tmp;
     max_val2 =(max_val2<fabs(tmp))?fabs(tmp):max_val2;
   }


      MPI_Barrier(Partition::getComm());
      MPI_Allreduce(&diff, &mean, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
      MPI_Allreduce(&diff2, &mean2, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
      MPI_Allreduce(&max_val, &max1, 1, MPI_DOUBLE, MPI_MAX, Partition::getComm());
      MPI_Allreduce(&max_val2,&max2, 1, MPI_DOUBLE, MPI_MAX, Partition::getComm());
      MPI_Allreduce(&num_halos, &n_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());
      MPI_Allreduce(&num_halos2, &n_tot2, 1, MPI_INT, MPI_SUM,  Partition::getComm());


   mean = mean/n_tot;
   mean2 = mean2/n_tot2;

   for (int i=0; i< num_halos; i++){
      double diff_tmp = (double)(val1->at(i))-mean;
      stddev += diff_tmp*diff_tmp/(n_tot-1);
   }
   for (int i=0; i< num_halos2; i++){
      double diff_tmp = (double)(val2->at(i))-mean2;
      stddev2 += diff_tmp*diff_tmp/(n_tot2-1);
   }


   MPI_Allreduce(&stddev, &stddev_tot, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   MPI_Allreduce(&stddev2, &stddev_tot2, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   stddev_tot = sqrt(stddev_tot);
   stddev_tot2 = sqrt(stddev_tot2);

   bool print_out = true;
   if (rank==0){

     if ((fabs((stddev_tot-stddev_tot2)/stddev_tot)<lim)&&(fabs((mean-mean2)/mean)<lim))
	 print_out=false;   
     if (isnan(mean)||isnan(mean2))
	 print_out=false;
     if (print_out){
     err ++;
     cout << var_name << endl;
     cout << var_name2 << endl;
     cout << "______________________________________" <<endl;
     cout << " means = " << mean << " and = "<< mean2 <<  endl;
     cout << " stddev = " << stddev_tot << " and  = " << stddev_tot2 << endl;
     cout << endl;
     }
   }
   return err;
}

int compute_mean_std_dist_halos2(Halos_test H_1 , Halos_test H_2, float lim ){
  // compute the mean and std of the differences 
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = H_1.num_halos;
  int count2 = H_2.num_halos;
  int err = 0;

  for (int i =0; i<N_HALO_FLOATS; i++){
    string var_name = float_var_names_test[i];
    string var_name2 = float_var_names_test2[i];
    err += compute_mean_float_dist(H_1.float_data[i],H_2.float_data[i],count,count2,var_name,var_name2, lim);
  }

  return err;
}


int compare_dist(string fof_file,string fof_file2, float lim){
	
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  if (rank==0){
     cout << " Comparing mean and standard deviations of property distributions" << endl;
     cout << " ________________________________________________________________" << endl;
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

  if (rank==0)
	  cout << endl;
  // Reading halos
  read_halos(H_1, fof_file, 2);
  read_halos(H_2, fof_file2, 2);

  if (rank == 0){
    cout << endl;
  }
  int n1 = H_1.num_halos;
  int n2 = H_2.num_halos;
  int n_tot, n_tot2;

  MPI_Allreduce(&n1, &n_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());
  MPI_Allreduce(&n2, &n_tot2, 1, MPI_INT, MPI_SUM,  Partition::getComm());


  // compute the error characteristics for the catalogs
  int err = compute_mean_std_dist_halos2(H_1 , H_2, lim);


    if ((rank==0)&&(err==0)){
      cout << " Results " << endl;
      cout << " _______ " << endl;
      cout << endl;
      cout << " Comparison test passed! " << endl;
      cout << " All variables within threshold of "  << lim << endl;
      cout << " Difference in number of halos  = "<< abs(n_tot-n_tot2) << endl;
      cout << endl;
  }
  if ((rank==0)&&(err>0)){
      cout << " Results " << endl;
      cout << " _______ " << endl;
      cout << endl;
      cout << " Comparison exceeded threshold of " << lim << " for " << err << " variables" << endl;
      cout << " out of a total of " <<  N_HALO_FLOATS << " variables " << endl;
      cout << " See above outputs for details  "<< endl;
      cout << " Difference in number of halos  = "<< abs(n_tot-n_tot2) << endl;
      cout << endl;
  }



  MPI_Barrier(Partition::getComm());
  H_1.Deallocate();
  H_2.Deallocate();


  return 0;
}
