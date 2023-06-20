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

#include <iostream>
#include <fstream>

#include <vector>
#include <climits>
#include <algorithm>

#include "H5Cpp.h"

#include "GenericIO.h"
#include "Halos_LJ.h"

// healpix header files
#include "chealpix.h"
#include "pointing.h"
#include "healpix_base.h"
#include "vec3.h"

#include "utils_pixel_LJ.h"

using namespace gio;
using namespace std;
using namespace H5;

#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


int main(int argc, char ** argv) {

    MPI_Init(&argc,&argv);
    int root_process = 0;
    int rank, numranks;
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 

    // input arguments
    char filename[512];
    strcpy(filename,argv[1]);
    long nside;
    nside = atoi(argv[2]);
    char zrange[5];
    strcpy(zrange,argv[3]);
    string zrange_str = zrange;

    // this is the filename base
    string base_name = filename;
    printf("filename = %s \n",filename);
    printf("basename = %s \n", base_name.c_str());


    // healpix setup
    Healpix_Ordering_Scheme ring = RING;
    Healpix_Ordering_Scheme nest = NEST;
    T_Healpix_Base<int> map_lores;  
    T_Healpix_Base<int> map_hires;
    int nside_hi = 1024;
    map_lores = Healpix_Base (nside, ring, SET_NSIDE);
    map_hires = Healpix_Base (nside_hi, ring, SET_NSIDE);
    long npix1 = map_lores.Npix();
    long npix_hi = map_hires.Npix();
    printf("number of pixels on sky = %ld \n", npix1);
    printf("number of pixels on sky (hires) = %ld \n", npix_hi);

    int status = 0;


    // work out which pixels are on the rank and initialize their output
    vector<int> pix_nums_rank;
    for (int j=0;j<npix1;++j){
       int rank_pix = compute_rank_fullsky(j, numranks);
       if (rank_pix==rank){
           pix_nums_rank.push_back(j);
       }
    }

      // rename this str_pix eventually 
      stringstream ss;
      ss<< rank;
      string str_pix = ss.str();
      cout << "/eagle/LastJourney/prlarsen/halo_lightcone_LJ/pixels_test4/" + zrange_str + "/cutout_"+str_pix+".hdf5"<< "\n";
      string new_name = "/eagle/LastJourney/prlarsen/halo_lightcone_LJ/pixels_test4/" + zrange_str + "/cutout_" + str_pix + ".hdf5";

      H5File* file = new H5File(new_name.c_str(), H5F_ACC_TRUNC);
      delete file;



    /*for (int l=0; l<pix_nums_rank.size(); l++){
      stringstream ss;
      ss<< pix_nums_rank[l];
      string str_pix = ss.str();
      cout << "/eagle/LastJourney/prlarsen/halo_lightcone_LJ/pixels_test2/" + zrange_str + "/cutout_"+str_pix+".hdf5"<< "\n";
      string new_name = "/eagle/LastJourney/prlarsen/halo_lightcone_LJ/pixels_test2/" + zrange_str + "/cutout_" + str_pix + ".hdf5";

      // create file and close  
      H5File* file = new H5File(new_name.c_str(), H5F_ACC_TRUNC);
      delete file;

    }*/



    // work out which steps were requested
    int num_files;
    vector<int> steps;
    int stepin;
    string filename_steps;
    filename_steps = "steps_total_"+zrange_str+".txt";
    ifstream myfile (filename_steps.c_str());
    if (myfile.is_open())
    {
      while (myfile >> stepin){
      steps.push_back(stepin);
      cout << stepin << '\n';
      }
      myfile.close();
    }

    
    for (int ifiles=0;ifiles<steps.size();ifiles++){
      
      Halos_test H_1; 
      H_1.Allocate();
      H_1.Set_MPIType();


      printf("ifiles = %d \n", ifiles);
      stringstream ss2;
      ss2<< steps[ifiles];
      string step_s = ss2.str();
      printf("step = %s \n", step_s.c_str());
      string mid = "/lc_intrp_halos_matched.";
      //string mid = "/lc_intrp_halos_fof_sod_matched.";
      string name = base_name + step_s + mid + step_s;
      printf("filename = %s \n", name.c_str());

      status = read_and_redistribute(name,map_lores,numranks, &H_1);

      vector<int64_t> pixel_counts;

      int pix_id_last = pix_nums_rank[0];
      int count_pix = 0;
      for (int j=0; j<H_1.num_halos;j++){
         int pix_id_tmp = H_1.pix_index->at(j);
         if (pix_id_tmp != pix_id_last){
            pixel_counts.push_back(j);
            count_pix++; 
            pix_id_last = pix_nums_rank[count_pix];
          }
      }
      pixel_counts.push_back(H_1.num_halos);

      assert(count_pix==pix_nums_rank.size()-1);
      //int64_t new_count=0;
      //for (int j=0;j<pix_nums_rank.size();j++){
      //    new_count += pixel_counts
      //}
      
      status = write_files_slab(H_1, rank, pix_nums_rank, pixel_counts, zrange_str, step_s);

      ///for (int l = 0; l< pix_nums_rank.size(); l++){
      //   status = write_files(H_1,  l, pix_nums_rank, pixel_counts, zrange_str, step_s);
      //}

    H_1.Deallocate();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

