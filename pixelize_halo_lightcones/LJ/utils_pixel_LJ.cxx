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
#include "chealpix.h"
#include "pointing.h"
#include "healpix_base.h"
#include "vec3.h"

#include "H5Cpp.h"
#include "H5Attribute.h"

#include "Halos_LJ.h"
#include "GenericIO.h"

#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


using namespace std;
using namespace gio;


//halo_properties_test
bool comp_rank(const halo_properties &a, const halo_properties &b){
   return a.rank < b.rank;
}

bool comp_pix(const halo_properties &a, const halo_properties &b){
   return a.pix_index < b.pix_index;
}


int compute_rank_fullsky(int pix_val, int numranks){
    int rankn = pix_val%numranks;
    return rankn;
}



void read_halos(Halos_test* H0, string file_name) {
  GenericIO GIO(MPI_COMM_WORLD,file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();

  (*H0).Resize(num_elems + GIO.requestedExtraSpace());

  GIO.addVariable("fof_halo_tag",   *((*H0).fof_halo_tag),true);

  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names[i], *((*H0).float_data[i]), true);
   
  for (int i=0; i<N_LIGHTCONE_FLOATS; ++i)
    GIO.addVariable((const string)lightcone_var_names[i], *((*H0).lightcone_data[i]), true);
 
  for (int i=0; i<N_HALO_FLOATS_E; ++i)
    GIO.addVariable((const string)float_var_names_ellipticity[i], *((*H0).ellipticity_data[i]), true);

  GIO.readData();
  (*H0).Resize(num_elems);
}


int redistribute_halos(int numranks, vector<halo_properties> halo_send, Halos_test* H_1, vector<int> send_count){

    vector<int> send_offset;
    send_offset.resize(numranks);
    send_offset[0] = 0;

    vector<int> recv_offset;
    recv_offset.resize(numranks);
    recv_offset[0] = 0;

    vector<int> recv_count;
    recv_count.resize(numranks);

    vector<halo_properties> recv_halos;

    MPI_Alltoall( &send_count[0], 1 , MPI_INT, &recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);

    for (int k=1;k<numranks;k++){
        send_offset[k] = send_offset[k-1] + send_count[k-1];
        recv_offset[k] = recv_offset[k-1] + recv_count[k-1];
    }
    recv_halos.resize(recv_offset.back()+recv_count.back());

    MPI_Alltoallv(&halo_send[0], &send_count[0], &send_offset[0],(*H_1).halo_properties_MPI_Type,\
                  &recv_halos[0], &recv_count[0], &recv_offset[0], (*H_1).halo_properties_MPI_Type, MPI_COMM_WORLD );

    int64_t nparticles = recv_offset.back()+recv_count.back();
    printf("nparticles = %d\n",nparticles);

    sort(recv_halos.begin(),recv_halos.end(),comp_pix);

    (*H_1).Resize(nparticles);
    for (int k=0; k<nparticles; k++){
        (*H_1).Assign(recv_halos[k],k);
    }
   

    return 0;
}


int compute_ranks_halo( Halos_test* H_1, T_Healpix_Base<int> map_lores, int numranks, vector<int> &send_count, vector<halo_properties>* halo_send){
 
     for (int k=0;k<(*H_1).num_halos;k++){
       int rankn;

       halo_properties tmp = (*H_1).GetProperties(k);
       // lightcone_data 0,1,2 are x,y,z
       vec3 vec_temp3 = vec3(tmp.lightcone_data[0],tmp.lightcone_data[1],tmp.lightcone_data[2]);
       int ind2 = map_lores.vec2pix(vec_temp3);
       rankn = compute_rank_fullsky(ind2, numranks);
      
       tmp.pix_index = ind2; 
       tmp.rank = rankn; 

       send_count[rankn] += 1;
       (*halo_send).push_back(tmp);
     
   }
 return 0;
}



int read_and_redistribute(string file_name, T_Healpix_Base<int> map_lores, int numranks, Halos_test* H_1){

    int status; 
    vector<halo_properties> halo_send;
    vector<int> send_count(numranks,0);


    read_halos(H_1, file_name); // read halos into H0 
    status = compute_ranks_halo(H_1, map_lores, numranks, send_count, &halo_send);
    sort(halo_send.begin(),halo_send.end(),comp_rank);
    status = redistribute_halos(numranks, halo_send, H_1, send_count);

    return status; 

}

int write_files_slab(Halos_test &H_1, int rank,  vector<int> pixel_nums_rank, vector<int64_t> pixel_counts, string zrange_str, string step_s){


        stringstream ss;
        ss<<rank;
        string str_rank = ss.str();

        string new_name = "/eagle/LastJourney/prlarsen/halo_lightcone_LJ/pixel_outputs/" + zrange_str + "/cutout_" + str_rank + ".hdf5";
        H5::H5File file(new_name.c_str(), H5F_ACC_RDWR); // must already exist 

        int64_t fill_value = 0;
        hsize_t chunk_dims[1] = {5000};
        H5::DSetCreatPropList plist;
        plist.setFillValue(H5::PredType::NATIVE_INT64, &fill_value);
        plist.setChunk(1,chunk_dims);

        float fill_value_float = 0;
        H5::DSetCreatPropList plist_float;
        plist_float.setFillValue(H5::PredType::NATIVE_FLOAT, &fill_value_float);
        plist_float.setChunk(1,chunk_dims);

        hsize_t dim_tot[1];
        dim_tot[0] = H_1.num_halos;
        H5::DataSpace dataspace(1,dim_tot); // number of dimensions, size of array 
        H5::Group group = file.createGroup(step_s.c_str());

        string var_name = "fof_halo_tag"; 
        string path_name = "/" + step_s + "/" + var_name;
        H5::DataSet dataset = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_INT64,dataspace, plist);
 
        vector<int64_t> start_vec(pixel_nums_rank.size(),0); 
        vector<int64_t> count_vec(pixel_nums_rank.size(),0);
        

        for (int l=0;l<pixel_counts.size();l++){
        int64_t start_val = (l==0)?0:pixel_counts[l-1];
        int64_t end_val = pixel_counts[l];

        hsize_t dim[1];
        dim[0] = end_val - start_val;
        cout << pixel_nums_rank[l] << " for a size of " << dim[0] <<  "start val = " << start_val<< " and end val = " << end_val << endl;

        hsize_t start[1];
        hsize_t stride[1];
        hsize_t count[1];
        hsize_t block[1];

        start[0] = start_val;
        stride[0] = 1;
        count[0] = end_val - start_val; 
        block[0] = 1; 
        start_vec[l] = start[0];
        count_vec[l] = count[0];

        H5::DataSpace dataspace_slab(1,dim);
        dataspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);
        dataset.write(&H_1.fof_halo_tag->at(start_val),H5::PredType::NATIVE_INT64,dataspace_slab,dataspace);
        dataspace_slab.close();

        }

        for (int i=0; i<N_HALO_FLOATS;i++){
          var_name = float_var_names[i];
          path_name = "/" + step_s + "/" + var_name;
          H5::DataSpace dataspace_var(1,dim_tot); 
          H5::DataSet dataset_var = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_FLOAT,dataspace_var, plist_float);
          for (int l=0;l<pixel_counts.size();l++){
            hsize_t dim[1];
            hsize_t start[1];
            hsize_t stride[1];
            hsize_t count[1];
            hsize_t block[1];
 
            dim[0] = count_vec[l];
            start[0] = start_vec[l];
            stride[0] = 1;
            block[0] = 1; 
            count[0] = count_vec[l];
            int64_t start_val = start_vec[l];
           
            H5::DataSpace dataspace_slab(1,dim);
            dataspace_var.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);
            dataset_var.write(&H_1.float_data[i]->at(start_val),H5::PredType::NATIVE_FLOAT,dataspace_slab,dataspace_var);
            dataspace_slab.close();
          }
          dataspace_var.close();
          dataset_var.close();
        }
        for (int i=0; i<N_LIGHTCONE_FLOATS;i++){
          var_name = lightcone_var_names[i];
          path_name = "/" + step_s + "/" + var_name;
          H5::DataSpace dataspace_var(1,dim_tot);
          H5::DataSet dataset_var = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_FLOAT,dataspace_var, plist_float);
          for (int l=0;l<pixel_counts.size();l++){
            hsize_t dim[1];
            hsize_t start[1];
            hsize_t stride[1];
            hsize_t count[1];
            hsize_t block[1];

            dim[0] = count_vec[l];
            start[0] = start_vec[l];
            stride[0] = 1;
            block[0] = 1;
            count[0] = count_vec[l];
            int64_t start_val = start_vec[l];

            H5::DataSpace dataspace_slab(1,dim);
            dataspace_var.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);
            dataset_var.write(&H_1.lightcone_data[i]->at(start_val),H5::PredType::NATIVE_FLOAT,dataspace_slab,dataspace_var);
            dataspace_slab.close();
          }
          dataspace_var.close();
          dataset_var.close();
        }
        for (int i=0; i<N_HALO_FLOATS_E;i++){
          var_name = float_var_names_ellipticity[i];
          path_name = "/" + step_s + "/" + var_name;
          H5::DataSpace dataspace_var(1,dim_tot);
          H5::DataSet dataset_var = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_FLOAT,dataspace_var, plist_float);
          for (int l=0;l<pixel_counts.size();l++){
            hsize_t dim[1];
            hsize_t start[1];
            hsize_t stride[1];
            hsize_t count[1];
            hsize_t block[1];

            dim[0] = count_vec[l];
            start[0] = start_vec[l];
            stride[0] = 1;
            block[0] = 1;
            count[0] = count_vec[l];
            int64_t start_val = start_vec[l];

            H5::DataSpace dataspace_slab(1,dim);
            dataspace_var.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);
            dataset_var.write(&H_1.ellipticity_data[i]->at(start_val),H5::PredType::NATIVE_FLOAT,dataspace_slab,dataspace_var);
            dataspace_slab.close();
          }
          dataspace_var.close();
          dataset_var.close();
        }



        dataspace.close();
        dataset.close();

 
        hsize_t dim_att[1];
        dim_att[0] = (hsize_t)pixel_counts.size();
         

        path_name = "/" + step_s + "/pixel_id" ;
        H5::DataSpace dataspace_p(1,dim_att);
        H5::IntType datatype(H5::PredType::NATIVE_INT);
        datatype.setOrder( H5T_ORDER_LE );
        H5::DataSet dataset_p = file.createDataSet(path_name.c_str(),datatype,dataspace_p);
        dataset_p.write(&pixel_nums_rank[0],H5::PredType::NATIVE_INT);
        dataset_p.close();
        dataspace_p.close();

        path_name = "/" + step_s + "/start_idx" ;
        H5::DataSpace dataspace_s(1,dim_att); 
        H5::IntType datatype_s(H5::PredType::NATIVE_INT64);
        datatype_s.setOrder( H5T_ORDER_LE );
        H5::DataSet dataset_s = file.createDataSet(path_name.c_str(),datatype_s,dataspace_s);
        dataset_s.write(&start_vec[0],H5::PredType::NATIVE_INT64);
        dataset_s.close();
        dataspace_s.close();

        path_name = "/" + step_s + "/count_idx" ;
        H5::DataSpace dataspace_c(1,dim_att);
        H5::DataSet dataset_c = file.createDataSet(path_name.c_str(),datatype_s,dataspace_c);
        dataset_c.write(&count_vec[0],H5::PredType::NATIVE_INT64);
        dataset_c.close();
        dataspace_c.close();


        

        
return 0;

}



int write_files(Halos_test &H_1, int l, vector<int> pixel_nums_rank, vector<int64_t> pixel_counts, string zrange_str, string step_s){

        stringstream ss;
        ss<< pixel_nums_rank[l];
        string str_pix = ss.str();

        string new_name = "/eagle/LastJourney/prlarsen/halo_lightcone_LJ/pixels_test2/" + zrange_str + "/cutout_" + str_pix + ".hdf5";
        H5::H5File file(new_name.c_str(), H5F_ACC_RDWR); // must already exist 

        int64_t start_val = (l==0)?0:pixel_counts[l-1];
        int64_t end_val = pixel_counts[l];


        //cout << start_val << end_val << std::endl;
        hsize_t dim[1];
        dim[0] = end_val - start_val;
        cout << pixel_nums_rank[l] << " for a size of " << dim[0] <<  "start val = " << start_val<< " and end val = " << end_val << endl;

        // do something for fof_halo_tag 
        string var_name = "fof_halo_tag";
        string path_name = "/" + step_s + "/" + var_name; 
        // create dataspace
        H5::DataSpace dataspace(1,dim); // number of dimensions, size of array 
        // create dataset
        //string group_name = step_s;a
        H5::Group group = file.createGroup(step_s.c_str());
        H5::DataSet dataset = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_INT64,dataspace);

        //vector<int64_t> fof_tag(dim[0],0);
        //for (int i=0;i<dim[0];i++){
        //  fof_tag[i] = H_1.fof_halo_tag->at(start_val+i);
        //}
        //H5::DataSet dataset = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_INT64,dataspace);
        // write dataset
        //
        //auto start = H_1.fof_halo_tag.
        //auto first = (*H_1.fof_halo_tag).begin()+start_val;
        //auto last = (*H_1.fof_halo_tag).begin()+end_val;
        //vector<int64_t> fof_tags(first,last);
        //vector<int64_t> fof_tags(H_1.fof_halo_tag->begin()+start_val, H_1.fof_halo_tag->begin()+end_val);
        //dataset.write(&fof_tags,H5::PredType::NATIVE_INT64);
        dataset.write(&H_1.fof_halo_tag->at(start_val),H5::PredType::NATIVE_INT64);
        //dataset.write(&H_1.fof_halo_tag.at(start_val),H5::PredType::NATIVE_INT64);
        //dataset.write(&fof_tag,H5::PredType::NATIVE_INT64);

        
        for (int i=0; i<N_HALO_FLOATS;i++){
          var_name = float_var_names[i];
          path_name = "/" + step_s + "/" + var_name;
          H5::DataSpace dataspace(1,dim); 
          H5::DataSet dataset = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_FLOAT,dataspace);
          dataset.write(&H_1.float_data[i]->at(start_val),H5::PredType::NATIVE_FLOAT);

        }
        for (int i=0; i<N_LIGHTCONE_FLOATS; i++){
          var_name = lightcone_var_names[i];
          path_name = "/" + step_s + "/" + var_name;
          H5::DataSpace dataspace(1,dim); // number of dimensions, size of array 
          H5::DataSet dataset = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_FLOAT,dataspace);
          dataset.write(&H_1.lightcone_data[i]->at(start_val),H5::PredType::NATIVE_FLOAT);

        }
        for (int i=0; i<N_HALO_FLOATS_E; i++){
          var_name = float_var_names_ellipticity[i];
          path_name = "/" + step_s + "/" + var_name;
          H5::DataSpace dataspace(1,dim); // number of dimensions, size of array 
          H5::DataSet dataset = file.createDataSet(path_name.c_str(),H5::PredType::NATIVE_FLOAT,dataspace);
          dataset.write(&H_1.ellipticity_data[i]->at(start_val),H5::PredType::NATIVE_FLOAT);
        }
        
             

return 0;

}

