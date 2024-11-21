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
#include <sys/stat.h>
#include <cmath>

#include "map_def.h"
#include "GenericIO.h"


// Cosmotools
using namespace std;
using namespace gio;

int check_file(string filename){
    int commRank, commRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commRanks);
    int valid = 0;
    if((commRank == 0) && (access(filename.c_str(), R_OK) == 0)){//using short-circuit
      valid = 1;
    }
    MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD); // step is valid
    return valid;
}

size_t read_gio_header_padded(string file_name){
  GenericIO GIO(MPI_COMM_WORLD,file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems() + GIO.requestedExtraSpace()/sizeof(double);
  return num_elems;
}
size_t read_gio_header(string file_name){
  GenericIO GIO(MPI_COMM_WORLD,file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems(); //+ GIO.requestedExtraSpace()/sizeof(double);
  return num_elems;
}


void read_map_file_hydro(MapDataHydro* M, string file_name) {
  // reads in maps
  GenericIO GIO(MPI_COMM_WORLD,file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();

  // parallel read - dereferencing pointer to allocate extra space 
  (*M).Resize(num_elems + GIO.requestedExtraSpace()/sizeof(double));

  //(*M) = value , (*M).double_data[i] is a pointer to the data, *((*M).double_data[i]) should pass you the vector directly
  for (int i=0; i<N_MAPS_HYDRO; i++){
    GIO.addVariable((const string)map_names_hydro[i], *((*M).double_data[i]), true);
  }
  GIO.addVariable("idx", (*M).pix_index,true);
  GIO.readData();
  (*M).Resize(num_elems);
}

void read_map_file_go(MapData* M, string file_name) {
  // reads in maps
  GenericIO GIO(MPI_COMM_WORLD,file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  // parallel read
  (*M).Resize(num_elems + GIO.requestedExtraSpace()/sizeof(double));
  for (int i=0; i<N_MAPS; i++)
    GIO.addVariable((const string)map_names_go[i], *((*M).double_data[i]), true);
  GIO.addVariable("idx", (*M).pix_index,true);
  GIO.readData();
  (*M).Resize(num_elems);
}

void write_map_file(MapData *M, string file_name, string old_file_name){
  // create GIO scope
  GenericIO GIO(MPI_COMM_WORLD,old_file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  double PhysOrigin[3];
  double PhysScale[3];
  GIO.readPhysOrigin(PhysOrigin);
  GIO.readPhysScale(PhysScale);

  GenericIO NewGIO(MPI_COMM_WORLD,file_name);
  int64_t size_rank = (*M).npixels;

  NewGIO.setNumElems(size_rank);
  for (int d=0; d < 3; ++d){
          NewGIO.setPhysOrigin(PhysOrigin[d], d);
          NewGIO.setPhysScale(PhysScale[d], d);
  }

  (*M).Resize(num_elems + GIO.requestedExtraSpace()/sizeof(double));
  for (int i = 0; i < N_MAPS; ++i){
      NewGIO.addVariable((const string)map_names_go[i], *((*M).double_data[i]), true); // no positional info
  }
  NewGIO.addVariable("idx", (*M).pix_index, true);
  NewGIO.write();
  MPI_Barrier(MPI_COMM_WORLD);
}

void write_map_file_hydro(MapDataHydro *M, string file_name, string old_file_name){
  // create GIO scope
  GenericIO GIO(MPI_COMM_WORLD,old_file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  double PhysOrigin[3];
  double PhysScale[3];
  GIO.readPhysOrigin(PhysOrigin);
  GIO.readPhysScale(PhysScale);

  GenericIO NewGIO(MPI_COMM_WORLD,file_name);
  int64_t size_rank = (*M).npixels;

  NewGIO.setNumElems(size_rank);
  for (int d=0; d < 3; ++d){
          NewGIO.setPhysOrigin(PhysOrigin[d], d);
          NewGIO.setPhysScale(PhysScale[d], d);
  }

  (*M).Resize(num_elems + GIO.requestedExtraSpace()/sizeof(double));
  for (int i = 0; i < N_MAPS_HYDRO; ++i){
      NewGIO.addVariable((const string)map_names_hydro[i], *((*M).double_data[i]), true); // no positional info
  }
  NewGIO.addVariable("idx", (*M).pix_index, true);
  NewGIO.write();
  MPI_Barrier(MPI_COMM_WORLD);
}



int main( int argc, char** argv ) {

  MPI_Init(&argc, &argv);
  int commRank, commRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

  if(argc != 9) {
     if (commRank==0){
     fprintf(stderr,"USAGE: %s <inputpath> <inputfile> <folders> <hydro> <nside> <step_start> <step_end> <outfile> \n", argv[0]);
     }
     exit(-1);
  }

  char filepath[512];
  char filename[512];
  char outfile[512];
  strcpy(filepath, argv[1]);
  strcpy(filename, argv[2]);
  strcpy(outfile, argv[8]);

  const char *mpiioName_base = filepath;
  const char *mpiioName_file = filename;

  bool folders = atoi(argv[3]); // should be 0 if false
  bool hydro = atoi(argv[4]); // should be 0 if false
  int64_t nside = atoi(argv[5]);
  int step_start = atoi(argv[6]);
  int step_end = atoi(argv[7]);

  //
  int64_t nside_long = (int64_t)nside;
  int64_t expected_size = nside_long * nside_long * 12;

  string start_file;
  //CHECK 1:  start with a check that all the steps ran and gave an output
  for (int jj=step_start;jj<step_end+1;jj++){
    MPI_Barrier(MPI_COMM_WORLD);
    char step[10*sizeof(char)];
    char mpiioName[512];
    sprintf(step,"%d",jj);
    if (folders){
      strcpy(mpiioName, mpiioName_base); 
      strcat(mpiioName, "/step_");
      strcat(mpiioName, step);
      strcat(mpiioName, "/");
      strcat(mpiioName, mpiioName_file);
      strcat(mpiioName, step);
      strcat(mpiioName, ".gio");
    }
    else{
      strcpy(mpiioName, mpiioName_base);
      strcat(mpiioName, "/");
      strcat(mpiioName, mpiioName_file);
      strcat(mpiioName, step);
      strcat(mpiioName, ".gio");
    }
    if (check_file(mpiioName)==0){
      if (commRank==0){
        fprintf(stderr,"File does not exist for step %s \n", step);
	fprintf(stderr,"output of check file %d \n", check_file(mpiioName));
	fprintf(stderr,"Filename is %s \n", mpiioName);
      }
      exit(-1);
    }
    assert(check_file(mpiioName)); // assert that the check worked
    if (jj==step_start){
      start_file = mpiioName;
    }
  } // close step loop

  if (commRank==0){
    fprintf(stdout,"Passed file existence check for all steps \n");
  }

  // allocating vectors so these stay in scope for now
  MapDataHydro M, M_sum;
  MapData MG, MG_sum;
  int64_t nelems_padded = read_gio_header_padded(start_file);
  int64_t nelems = read_gio_header(start_file);

  if (hydro){
	  M.Allocate(nelems_padded);
	  M_sum.Allocate(nelems); // allocating extra space for GIO read, will resize back down after read
  }
  else{
	  MG.Allocate(nelems_padded);
	  MG_sum.Allocate(nelems);
  }
  
  ofstream outfile_min, outfile_max, outfile_std, outfile_mean, outfile_count0, outfile_countbad;
  if (commRank==0){
  outfile_min.open("logs/logfile_min.txt");
  outfile_max.open("logs/logfile_max.txt");
  outfile_std.open("logs/logfile_std.txt"); 
  outfile_mean.open("logs/logfile_mean.txt"); 
  outfile_count0.open("logs/logfile_count0.txt"); 
  outfile_countbad.open("logs/logfile_countbad.txt"); 
  }

  

  //TODO: assert that the pixel indices are the same for each step. Initialize that at the first step then add a check for every subsequent step
 
  // CHECKS 2-5 require a map read, confirming number of pixels, outputting mean, stddev, min, max, nzeros, nbad, creating summed map
  for (int jj=step_start;jj<step_end+1;jj++){
    MPI_Barrier(MPI_COMM_WORLD);
    char step[10*sizeof(char)];
    char mpiioName[512];
    sprintf(step,"%d",jj);
    if (folders){
      strcpy(mpiioName, mpiioName_base);
      strcat(mpiioName, "/step_");
      strcat(mpiioName, step);
      strcat(mpiioName, "/");
      strcat(mpiioName, mpiioName_file);
      strcat(mpiioName, step);
      strcat(mpiioName, ".gio");
    }
    else{
      strcpy(mpiioName, mpiioName_base);
      strcat(mpiioName, "/");
      strcat(mpiioName, mpiioName_file);
      strcat(mpiioName, step);
      strcat(mpiioName, ".gio");
    }
    assert(check_file(mpiioName)); // assert that the file exists

    if (hydro){
      // check on number of pixels in the file 
      read_map_file_hydro(&M, mpiioName);
      int64_t size_global;
      int64_t size_rank = M.npixels;
      MPI_Reduce(&size_rank,&size_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&size_global,1,MPI_INT64_T, 0, MPI_COMM_WORLD);

      if (size_global!=expected_size){
        if (commRank==0){
          fprintf(stderr,"For step %s, total pixel count is %ld and expected count is %ld  \n", step, size_global, expected_size); // not sure if I need a long here
        }
        M.Deallocate();
        M_sum.Deallocate();
        MPI_Finalize();
        exit(-1);
      }
      assert(size_global==expected_size);

      if (jj==step_start){
        M_sum.Clear(size_rank);
        for (size_t kk=0;kk<M_sum.npixels;kk++){
          M_sum.pix_index.at(kk) = M.pix_index.at(kk);
        }
      }

      for (int ii=0;ii<N_MAPS_HYDRO;ii++){
        double min_val_global;
        double max_val_global;
        double sum_val_global;
        double std_val_global;
        double mean_val_global;
        double min_val = 1.e10;
        double max_val = -1.e10;
        double sum_val = 0.0;
        int64_t count_zero = 0;
        int64_t count_bad = 0;
        int64_t count_zero_global;
        int64_t count_bad_global;
	int64_t bad_pixels=0;
	int64_t bad_pixels_global;

        for (int64_t pp=0;pp<M.npixels;pp++){
	  double val = M.double_data[ii]->at(pp);
	  min_val = (val<min_val)?val:min_val;
	  max_val = (val>max_val)?val:max_val;
	  sum_val += val;
	  if (val==0){
	    count_zero++;
	  } // assume we have cmath
	  if (isnan(val) or isinf(val)){
            count_bad++;
	  }
	  M_sum.double_data[ii]->at(pp) += M.double_data[ii]->at(pp);
          if (M_sum.pix_index.at(pp)!=M.pix_index.at(pp)){
              bad_pixels +=1;
          }
        }
        MPI_Reduce(&bad_pixels,&bad_pixels_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Bcast(&bad_pixels_global,1,MPI_INT64_T, 0, MPI_COMM_WORLD);
        if (bad_pixels_global>0){
          if (commRank==0){
            fprintf(stderr,"For step %s, pixel mismatch relative to starting step for %ld pixels  \n", step, bad_pixels_global); // not sure if I need a long here
          }
          M.Deallocate();
          M_sum.Deallocate();
          MPI_Finalize();
          exit(-1);
        }


        MPI_Reduce(&max_val,&max_val_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        MPI_Reduce(&sum_val,&sum_val_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&min_val,&min_val_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
        MPI_Reduce(&count_zero,&count_zero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&count_bad,&count_bad_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);

        if (commRank==0){
          mean_val_global = sum_val_global/(double)size_global;
        }
        MPI_Bcast(&mean_val_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // compute std dev 
        sum_val = 0.0;
        for (int64_t pp=0;pp<M.npixels;pp++){
          double val = M.double_data[ii]->at(pp);
	  sum_val += (val-mean_val_global)*(val-mean_val_global);
        }
        MPI_Reduce(&sum_val,&std_val_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        std_val_global = std_val_global/(double)size_global;
        std_val_global = sqrt(std_val_global); // assuming using namespace std
        // if commRank=0 put the values into a text file 

        if (commRank==0){
          if (ii<N_MAPS_HYDRO-1){
            outfile_min << min_val_global << ",";
            outfile_max << max_val_global << ",";
            outfile_std << std_val_global << ",";
            outfile_mean << mean_val_global << ",";
            outfile_count0 << count_zero_global << ",";
            outfile_countbad << count_bad_global << ",";
          }
          else{
            outfile_min << min_val_global << "\n";
            outfile_max << max_val_global << "\n";
            outfile_std << std_val_global << "\n";
            outfile_mean << mean_val_global << "\n";
            outfile_count0 << count_zero_global << "\n";
            outfile_countbad << count_bad_global << "\n";
          }
        }
      } // closing over N_MAPS_HYDRO
    } // closing over hydro 
    else{
      // check on number of pixels in the file 
      read_map_file_go(&MG, mpiioName);
      int64_t size_global;
      int64_t size_rank = MG.npixels;
      MPI_Reduce(&size_rank,&size_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&size_global,1,MPI_INT64_T, 0, MPI_COMM_WORLD);
      if (size_global!=expected_size){
        if (commRank==0){
          fprintf(stderr,"For step %s, total pixel count is %ld and expected count is %ld  \n", step, size_global, expected_size); // not sure if I need a long here
        }
	MG.Deallocate();
        MG_sum.Deallocate();
        MPI_Finalize();
        exit(-1);
      }
      assert(size_global==expected_size);
      MPI_Barrier(MPI_COMM_WORLD);

       if (commRank==0){
          fprintf(stdout,"For step %s, total pixel count is correct: %ld and expected count is %ld  \n", step, size_global, expected_size); // not sure if I need a long here
        }


      // initialize M_sum to 0s here 
      if (jj==step_start){
        MG_sum.Clear(size_rank);
        for (size_t kk=0;kk<MG_sum.npixels;kk++){
          MG_sum.pix_index.at(kk) = MG.pix_index.at(kk);
        }
      }

      for (int ii=0;ii<N_MAPS;ii++){
        double min_val_global;
        double max_val_global;
        double sum_val_global;
        double std_val_global;
        double mean_val_global;
        double min_val = 1.e10;
        double max_val = -1.e10;
        double sum_val = 0.0;
        int64_t count_zero = 0;
        int64_t count_bad = 0;
        int64_t count_zero_global;
        int64_t count_bad_global;
	int64_t bad_pixels = 0;
	int64_t bad_pixels_global;

        for (int64_t pp=0;pp<MG.npixels;pp++){
          double val = MG.double_data[ii]->at(pp);
          min_val = (val<min_val)?val:min_val;
          max_val = (val>max_val)?val:max_val;
          sum_val += val;
          if (val==0){
            count_zero++;
          } // assume we have cmath
          if (isnan(val) or isinf(val)){
            count_bad++;
          }
          MG_sum.double_data[ii]->at(pp) += MG.double_data[ii]->at(pp);
          if (MG_sum.pix_index.at(pp)!=MG.pix_index.at(pp)){
              bad_pixels +=1;
	  }
        }
        MPI_Reduce(&bad_pixels,&bad_pixels_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Bcast(&bad_pixels_global,1,MPI_INT64_T, 0, MPI_COMM_WORLD);
        if (bad_pixels_global>0){
          if (commRank==0){
            fprintf(stderr,"For step %s, pixel mismatch relative to starting step for %ld pixels  \n", step, bad_pixels_global); // not sure if I need a long here
          }
          MG.Deallocate();
          MG_sum.Deallocate();
          MPI_Finalize();
          exit(-1);
        }



        MPI_Reduce(&max_val,&max_val_global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        MPI_Reduce(&sum_val,&sum_val_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&min_val,&min_val_global,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
        MPI_Reduce(&count_zero,&count_zero_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&count_bad,&count_bad_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);

        if (commRank==0){
          mean_val_global = sum_val_global/(double)size_global;
        }
        MPI_Bcast(&mean_val_global,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // compute std dev 
        sum_val = 0.0;
        for (int64_t pp=0;pp<MG.npixels;pp++){
          double val = MG.double_data[ii]->at(pp);
          sum_val += (val-mean_val_global)*(val-mean_val_global);
        }
        MPI_Reduce(&sum_val,&std_val_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        std_val_global = std_val_global/(double)size_global;
        std_val_global = sqrt(std_val_global); // assuming using namespace std
        // if commRank=0 put the values into a text file 

        if (commRank==0){
          if (ii<N_MAPS-1){
            outfile_min << min_val_global << ",";
            outfile_max << max_val_global << ",";
            outfile_std << std_val_global << ",";
            outfile_mean << mean_val_global << ",";
            outfile_count0 << count_zero_global << ",";
            outfile_countbad << count_bad_global << ",";
          }
          else{
            outfile_min << min_val_global << "\n";
            outfile_max << max_val_global << "\n";
            outfile_std << std_val_global << "\n";
            outfile_mean << mean_val_global << "\n";
            outfile_count0 << count_zero_global << "\n";
            outfile_countbad << count_bad_global << "\n";
          }
        }
      } // closes over N_MAPS


    } // closes over if/else for hydro 
  }// close over steps

  if (commRank==0){
    outfile_min.close(); // update to an input name
    outfile_max.close();
    outfile_std.close();
    outfile_mean.close();
    outfile_count0.close();
    outfile_countbad.close();
  }
  // write sum to a file 
  if (hydro){
    write_map_file_hydro(&M_sum, outfile, start_file);
  }
  else{
    write_map_file(&MG_sum, outfile, start_file);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (hydro){
  M.Deallocate();
  M_sum.Deallocate();
  }
  else{
  MG.Deallocate();
  MG_sum.Deallocate();
  }
  
  MPI_Finalize();
  return 0;
}
