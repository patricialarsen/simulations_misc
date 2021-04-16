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


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;


//vector<float> x;


int find_log_bin(float value, float min_val, float bin_width){
    int bin_idx = (int)((log(value)-log(min_val))/bin_width);
    return bin_idx;
}


int find_linear_bin(float value, float min_val, float bin_width){
    int bin_idx = (int)((value-min_val)/bin_width);
    return bin_idx;
}

void read_halo_file(string file_name, string var, vector<float> &x) {
  {
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  x.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  GIO.addVariable(var, x, true);
  GIO.readData();
  x.resize(num_elems);
  }
  return; 
}

void write_bins_2d(string file_name, vector<float> bin_count_tot, int nbins){
  ofstream file_bins (file_name);
  if (file_bins.is_open())
  {
   for (int i=0;i<nbins;i++){
       for (int j=0; j<nbins; j++){
       file_bins << bin_count_tot[i,j]<< " ";
      }
       file_bins <<"\n";
    }
    file_bins.close();
   }
  return;
}



void write_bins(string file_name, vector<float> bin_count_tot, int nbins){
  ofstream file_bins (file_name);
  if (file_bins.is_open())
  {
   for (int i=0;i<nbins;i++){
       file_bins << bin_count_tot[i] << "\n";
    }
    file_bins.close();
   }
  return;
}

void make_mass_hist(vector<float> x, vector<float> m, float m_min, float m_max, int nbins, float min_val, float max_val, string filename, int rank){
   int count = (int)x.size();
   vector<int> bin_count(nbins,nbins);
   
   float bin_width = (log(max_val)-log(min_val))/nbins;
   float bin_width_mass = (log(m_max)- log(m_min))/nbins;
   for (int i=0; i< count; ++i){
      int idx_tmp1 = find_log_bin(x[i],min_val, bin_width);
      int idx_tmp2 = find_log_bin(m[i],m_min, bin_width_mass);
      if ((idx_tmp1>=0)&&(idx_tmp1<nbins)&&(idx_tmp2>0)&&(idx_tmp2<nbins)){
          bin_count[idx_tmp1,idx_tmp2]+=1;
      }
   }


  int count_tot;
  vector<int> bin_count_tot(nbins,nbins);
  vector<float> bin_count_tot_float(nbins,nbins);

  MPI_Reduce(&count,&count_tot,1,MPI_INT,MPI_SUM,0,Partition::getComm());
  for (int i=0; i<nbins; i++){
  MPI_Reduce(&bin_count[i,0],&bin_count_tot[i,0],nbins,MPI_INT,MPI_SUM,0,Partition::getComm());
  }
  if (rank == 0){
     for (int i=0; i<nbins; i++){
        for (int j=0; j<nbins;j++){
        bin_count_tot_float[i,j] = (float)bin_count_tot[i,j]/count_tot;
       }
     }
     write_bins_2d(filename,bin_count_tot_float,nbins);
  }


}



void make_linear_histogram_mbins(vector<float> x, vector<float> m, float m_min, float m_max, int nbins, float min_val, float max_val, string filename, int rank){

  double sum_x=0;
  int count = (int)x.size();
  vector<int> bin_count(nbins,0);

  float bin_width = (max_val-min_val)/nbins;

  for (int i=0; i<count; ++i) {
      int idx_tmp = find_linear_bin(x[i], min_val,bin_width);
      if ((idx_tmp>=0)&&(idx_tmp<nbins)&&(m[i]<m_max)&&(m[i]>=m_min)){
        bin_count[idx_tmp]+=1;
      }
      sum_x += (double)x[i];
  }

  double sum_tot;
  int count_tot;
  vector<int> bin_count_tot(nbins);
  vector<float> bin_count_tot_float(nbins);

  MPI_Reduce(&sum_x,&sum_tot,1,MPI_DOUBLE,MPI_SUM,0,Partition::getComm());
  MPI_Reduce(&count,&count_tot,1,MPI_INT,MPI_SUM,0,Partition::getComm());
  MPI_Reduce(&bin_count[0],&bin_count_tot[0],nbins,MPI_INT,MPI_SUM,0,Partition::getComm());

  if (rank == 0){
     double mean_tot = sum_tot/count_tot;
     cout << mean_tot << endl;
     for (int i=0; i<nbins; i++){
        bin_count_tot_float[i] = (float)bin_count_tot[i]/count_tot;
      }
     write_bins(filename,bin_count_tot_float,nbins);
  }

}

void make_log_histogram(vector<float> x, int nbins, float min_val, float max_val, string filename, int rank){

  double sum_x=0;
  int count = (int)x.size();
  vector<int> bin_count(nbins,0);

  float bin_width = (log(max_val)-log(min_val))/nbins;

  for (int i=0; i<count; ++i) {
      int idx_tmp = find_log_bin(x[i], min_val,bin_width);
      if ((idx_tmp>=0)&&(idx_tmp<nbins)){
        bin_count[idx_tmp]+=1;
      }
      sum_x += (double)x[i];
  }

  double sum_tot;
  int count_tot;
  vector<int> bin_count_tot(nbins);
  vector<float> bin_count_tot_float(nbins);

  MPI_Reduce(&sum_x,&sum_tot,1,MPI_DOUBLE,MPI_SUM,0,Partition::getComm());
  MPI_Reduce(&count,&count_tot,1,MPI_INT,MPI_SUM,0,Partition::getComm());
  MPI_Reduce(&bin_count[0],&bin_count_tot[0],nbins,MPI_INT,MPI_SUM,0,Partition::getComm());

  if (rank == 0){
     double mean_tot = sum_tot/count_tot;
     cout << mean_tot << endl;
     for (int i=0; i<nbins; i++){
        bin_count_tot_float[i] = (float)bin_count_tot[i]/count_tot;
      }
     write_bins(filename,bin_count_tot_float,nbins);
  }


}




void make_linear_histogram(vector<float> x, int nbins, float min_val, float max_val, string filename, int rank){

  double sum_x=0;
  int count = (int)x.size();
  vector<int> bin_count(nbins,0);

  float bin_width = (max_val-min_val)/nbins;

  for (int i=0; i<count; ++i) {
      int idx_tmp = find_linear_bin(x[i], min_val,bin_width);
      if ((idx_tmp>=0)&&(idx_tmp<nbins)){
        bin_count[idx_tmp]+=1;
      }
      sum_x += (double)x[i];
  }

  double sum_tot;
  int count_tot;
  vector<int> bin_count_tot(nbins);
  vector<float> bin_count_tot_float(nbins);

  MPI_Reduce(&sum_x,&sum_tot,1,MPI_DOUBLE,MPI_SUM,0,Partition::getComm());
  MPI_Reduce(&count,&count_tot,1,MPI_INT,MPI_SUM,0,Partition::getComm());
  MPI_Reduce(&bin_count[0],&bin_count_tot[0],nbins,MPI_INT,MPI_SUM,0,Partition::getComm());

  if (rank == 0){
     double mean_tot = sum_tot/count_tot;
     cout << mean_tot << endl;
     for (int i=0; i<nbins; i++){
        bin_count_tot_float[i] = (float)bin_count_tot[i]/count_tot;
      }
     write_bins(filename,bin_count_tot_float,nbins);
  }


}

void make_spins(vector<float> x, vector<float> y, vector<float> z, vector<float> m, vector<float> &s){
    s.resize(x.size());
    for (int i;i<x.size();i++){
       s[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i])/m[i];
    }
}



void angularmomentum_histograms(string fof_file, float ammax, int nbins, int rank){
  vector<float> x;
  vector<float> y;
  vector<float> z;
  vector<float> m;
  read_halo_file(fof_file,"fof_halo_angmom_x", x);
  read_halo_file(fof_file,"fof_halo_angmom_y", y);
  read_halo_file(fof_file,"fof_halo_angmom_z", z);
  read_halo_file(fof_file,"fof_halo_mass", m);


  if (rank == 0)
    cout << "Read fof halos" << endl;
  MPI_Barrier(Partition::getComm());

  float min_val = -ammax;
  float max_val = ammax;
  float m_min = 1.e12;
  float m_max = 1.e15;
  make_linear_histogram_mbins(x, m, m_min, m_max, nbins, min_val, max_val, "outputs/amx_hist.txt", rank);
  make_linear_histogram_mbins(y, m, m_min, m_max,  nbins, min_val, max_val, "outputs/amy_hist.txt", rank);
  make_linear_histogram_mbins(z, m, m_min, m_max, nbins, min_val, max_val, "outputs/amz_hist.txt", rank);

  vector<float> s;
  make_spins(x,y,z,m,s);
  make_log_histogram(s,nbins,1.e-2,10.,"outputs/spin_hist.txt",rank);
  make_mass_hist(s,m,m_min,m_max, nbins, 1.e-2, 10., "outputs/spin_hist2d.txt",rank);
}


void velocity_histograms(string fof_file, float vmax, int nbins, int rank){
  vector<float> x;
  vector<float> y;
  vector<float> z;
  read_halo_file(fof_file,"fof_halo_com_vx", x);
  read_halo_file(fof_file,"fof_halo_com_vy", y);
  read_halo_file(fof_file,"fof_halo_com_vz", z);

  if (rank == 0)
    cout << "Read fof halos" << endl;
  MPI_Barrier(Partition::getComm());

  float min_val = -vmax;
  float max_val = vmax;
  make_linear_histogram(x, nbins, min_val, max_val, "outputs/vx_hist.txt", rank);
  make_linear_histogram(y, nbins, min_val, max_val, "outputs/vy_hist.txt", rank);
  make_linear_histogram(z, nbins, min_val, max_val, "outputs/vz_hist.txt", rank);

}


void position_histograms(string fof_file, float boxsize, int nbins, int rank){
  vector<float> x;
  vector<float> y;
  vector<float> z;
  read_halo_file(fof_file,"fof_halo_center_x", x);
  read_halo_file(fof_file,"fof_halo_center_y", y);
  read_halo_file(fof_file,"fof_halo_center_z", z);

  if (rank == 0)
    cout << "Read fof halos" << endl;
  MPI_Barrier(Partition::getComm());

  float min_val = 0;
  float max_val = boxsize;
  make_linear_histogram(x, nbins, min_val, max_val, "outputs/x_hist.txt", rank);
  make_linear_histogram(y, nbins, min_val, max_val, "outputs/y_hist.txt", rank);
  make_linear_histogram(z, nbins, min_val, max_val, "outputs/z_hist.txt", rank);

}

int main( int argc, char** argv ) {
  MPI_Init( &argc, &argv );
  Partition::initialize();
  GenericIO::setNaturalDefaultPartition();

  int rank, n_ranks;
  int nbins = 100;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  string fof_file    = string(argv[1]);

  position_histograms(fof_file, 256.,nbins,rank);
  velocity_histograms(fof_file,1000.,nbins,rank);
  angularmomentum_histograms(fof_file,1.e16,nbins,rank);

  Partition::finalize();
  MPI_Finalize();
  return 0;
}
