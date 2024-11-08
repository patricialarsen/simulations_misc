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

void read_map_file_hydro(MapDataHydro* M, string file_name) {
  // reads in maps
  GenericIO GIO(MPI_COMM_WORLD,file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  cout << "this is the number of elements on my rank " << num_elems << endl;
  // parallel read
  (*M).Resize(num_elems + GIO.requestedExtraSpace()/sizeof(double));
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
  (*M).Resize(num_elems + GIO.requestedExtraSpace());
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
  vector<vector<double>> var_data_doubles(N_MAPS);
  for (int i = 0; i < N_MAPS; ++i){
      var_data_doubles[i].resize(size_rank + NewGIO.requestedExtraSpace()/sizeof(double));
      for (int64_t j=0; j<size_rank; j++){
          var_data_doubles[i][j] =  (*M->double_data[i])[j];
      }
      NewGIO.addVariable((const string)map_names_go[i], var_data_doubles[i], true); // no positional info
  }
  (*M).pix_index.resize(size_rank + NewGIO.requestedExtraSpace()/sizeof(double));
  NewGIO.addVariable("idx", (*M).pix_index, true);
  NewGIO.write();
  (*M).pix_index.resize(size_rank); // shouldn't be strictly necessary but for safety
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
  vector<vector<double>> var_data_doubles(N_MAPS_HYDRO);
  for (int i = 0; i < N_MAPS; ++i){
      var_data_doubles[i].resize(size_rank + NewGIO.requestedExtraSpace()/sizeof(double));
      for (int64_t j=0; j<size_rank; j++){
          var_data_doubles[i][j] =  (*M->double_data[i])[j];
      }
      NewGIO.addVariable((const string)map_names_hydro[i], var_data_doubles[i], true); // no positional info
  }
  (*M).pix_index.resize(size_rank + NewGIO.requestedExtraSpace()/sizeof(double));
  NewGIO.addVariable("idx", (*M).pix_index, true);
  NewGIO.write();
  (*M).pix_index.resize(size_rank); // shouldn't be strictly necessary but for safety
  MPI_Barrier(MPI_COMM_WORLD);
}



int main( int argc, char** argv ) {

  MPI_Init(&argc, &argv);
  int commRank, commRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commRanks);

  if(argc != 8) {
     if (commRank==0){
     fprintf(stderr,"USAGE: %s <inputpath> <inputfile> <folders> <hydro> <nside> <step_start> <step_end> \n", argv[0]);
     }
     exit(-1);
  }

  //TODO: add output file path for summed data 
  char filepath[512];
  char filename[512];
  strcpy(filepath,argv[1]);
  strcpy(filename,argv[2]);
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
  } // close step loop

  if (commRank==0){
    fprintf(stdout,"Passed file existence check for all steps \n");
  }

  // allocating vectors so these stay in scope for now
  MapDataHydro M, M_sum;
  MapData MG, MG_sum;

  M.Allocate();
  MG.Allocate();
  M_sum.Allocate();
  MG_sum.Allocate();

  ofstream outfile_min, outfile_max, outfile_std, outfile_mean, outfile_count0, outfile_countbad;
  if (commRank==0){
  outfile_min.open("logfile_min.txt"); // update to an input name
  outfile_max.open("logfile_max.txt");
  outfile_std.open("logfile_std.txt"); 
  outfile_mean.open("logfile_mean.txt"); 
  outfile_count0.open("logfile_count0.txt"); 
  outfile_countbad.open("logfile_countbad.txt"); 
  }


  // assert that the pixel indices are the same for each step. Initialize that at the first step then add a check for every subsequent step
 
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
      if (size_global!=expected_size){
        if (commRank==0){
          fprintf(stderr,"For step %s, total pixel count is %ld and expected count is %ld  \n", step, size_global, expected_size); // not sure if I need a long here
        }
        exit(-1);
      }
      assert(size_global==expected_size);

      if (jj==step_start){
      M_sum.Clear(size_rank);
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
      int64_t size_rank = M.npixels;
      MPI_Reduce(&size_rank,&size_global,1,MPI_INT64_T,MPI_SUM,0,MPI_COMM_WORLD);
      if (size_global!=expected_size){
        if (commRank==0){
          fprintf(stderr,"For step %s, total pixel count is %ld and expected count is %ld  \n", step, size_global, expected_size); // not sure if I need a long here
        }
        exit(-1);
      }
      assert(size_global==expected_size);

      // initialize M_sum to 0s here 
      if (jj==step_start){
        MG_sum.Clear(size_rank);
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
          MG_sum.double_data[ii]->at(pp) += M.double_data[ii]->at(pp);
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

  MPI_Barrier(MPI_COMM_WORLD);
  M.Deallocate();
  MG.Deallocate();
  M_sum.Deallocate();
  MG_sum.Deallocate();
  MPI_Finalize();
  return 0;
}
