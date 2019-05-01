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

// dtk and pixelize header files
#include "pixelize_halos_sod.h"
#include "hdf5_util.hpp"
#include "GenericIO.h"

// healpix header files
#include "chealpix.h"
#include "pointing.h"
#include "healpix_base.h"
#include "vec3.h"

#include "utils_pixel.h"

using namespace gio;
using namespace std;


#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


int main(int argc, char ** argv) {

    MPI_Init(&argc,&argv);
    int root_process = 0;
    int rank, numranks;
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
    MPI_Datatype particles_mpi =  createparticles();

    // input arguments
    char filename[512];
    strcpy(filename,argv[1]);
    long nside;
    nside = atoi(argv[2]);
    int octant=1;

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
    vector<int> pix_list_oct;
    vector<int> pix_list_noct;


    int status = 0;
    // work out pixels on octant (with some overlap)
    if (octant==1){
      double edge_val=0.1;
      status = compute_oct_pixels(map_hires, map_lores, pix_list_oct, pix_list_noct, edge_val);
     }

    // determine pixel numbers within this rank (just to create files)
    vector<int> pix_nums_rank;
    for (int j=0;j<pix_list_oct.size();++j){
       int rank_pix = compute_rank_partsky(pix_list_oct[j], numranks,pix_list_oct);
       if (rank_pix==rank){
           pix_nums_rank.push_back(pix_list_oct[j]);
       } 
    }
    // delete and recreate all required files
    for (int l=0; l<pix_nums_rank.size(); l++){
      stringstream ss;
      ss<< pix_nums_rank[l];
      string str_pix = ss.str();
      cout << "/projects/DarkUniverse_esp/prlarsen/healpix_cutouts_dc2_32_5000run_test/" + zrange_str + "/cutout_" + str_pix + ".hdf5" << "\n";
      H5::H5File file("/projects/DarkUniverse_esp/prlarsen/healpix_cutouts_dc2_32_5000run_test/" + zrange_str + "/cutout_" + str_pix + ".hdf5", H5F_ACC_TRUNC);
    }


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
    // opening file loop
      vector<int> send_count(numranks,0);
      vector<particles> particle_data1;
      vector<particles> recv_particles;
  
      printf("ifiles = %d \n", ifiles);
      stringstream ss2;
      ss2<< steps[ifiles];
      string step_s = ss2.str();
      printf("step = %s \n", step_s.c_str());
      string mid = "/lc_intrp_halos_fof_sod_matched.";
      string name = base_name + step_s + mid + step_s;
      printf("filename = %s \n", name.c_str());

      // deal with missing steps
      if (step_s=="347"){
      printf("hit step 347 - file doesn't exist, reading in 338 instead\n");
      name = base_name + "338" + mid + "338";
      printf("filename = %s \n", name.c_str());
      }

      else if (step_s=="259"){
      printf("hit step 259 - file doesn't exist, reading in 253 instead\n");
      name = base_name + "253" + mid + "253";
      printf("filename = %s \n", name.c_str());
      }

      int nparticles;
      printf("filename = %s \n", name.c_str());
      status = read_and_redistribute(name,map_lores,numranks,octant,send_count,recv_particles,particles_mpi,nparticles,pix_list_oct, step_s);
   

      int i,j;
      // determine which pixels actually contain particles within this rank
      vector<int> pixel_ids;
      vector<int> pixel_counts;

      for (i=0;i<pix_list_oct.size();i++){
         int rank_pix = compute_rank_partsky(pix_list_oct[i], numranks, pix_list_oct);
         if (rank_pix==rank){
            pixel_ids.push_back(pix_list_oct[i]);
            int tmp = 0;
            for (j=0;j<nparticles;j++){
              int count_new;
              count_new = recv_particles[j].pix_index;
              if (count_new==pix_list_oct[i]){
                 tmp+=1;
              }
            }
            pixel_counts.push_back(tmp);
         }
      }

      for (int l = 0; l< pixel_counts.size(); l++){
          // write out hdf5 files for each of the pixels
          status = write_files(recv_particles, nparticles, l, pixel_ids, zrange_str, step_s);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

