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

#include "hdf5.h"
//#include "H5Cpp.h"
//#include "H5Attribute.h"

#include "MapDef.h"
#include "GenericIO.h"

#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


using namespace std;
using namespace gio;
//using namespace H5;

bool comp_rank(const map_properties &a, const map_properties &b){
     	return a.rank < b.rank;
}

bool comp_pix(const map_properties &a, const map_properties &b){
   return a.pix_index < b.pix_index;
}

/*bool comp_rank_hydro(const map_properties &a, const map_properties &b){
        return a.rank < b.rank;
}

bool comp_pix(const map_properties &a, const map_properties &b){
   return a.pix_index < b.pix_index;
}
*/


void read_map_file(MapData* M, string file_name) {
  // reads in maps
  GenericIO GIO(MPI_COMM_WORLD,file_name);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
    
  // parallel read - dereferencing pointer to allocate extra space 
  (*M).Resize(num_elems + GIO.requestedExtraSpace()/sizeof(double));
  
  //(*M) = value , (*M).double_data[i] is a pointer to the data, *((*M).double_data[i]) should pass you the vector directly
  for (int i=0; i<N_MAPS; i++){
    GIO.addVariable((const string)map_names_go[i], *((*M).double_data[i]), true);
  }   
  GIO.addVariable("idx", (*M).pix_index,true);
  GIO.readData();
  (*M).Resize(num_elems);
} 
  
void rescale_map_file(MapData* M){
    size_t npix = (*M).npixels;
    //TODO: hard-coded for Frontier-E simulation currently
    double MPP = 1342777241.6; 
    size_t npix_tot = (*M).npixels_tot;
    double pixsize = (4.*3.141529/npix_tot);
    // assert that the order is phi, vel, rho
    for (int64_t i=0;i<npix;i++){
      if ((*M).double_data[2]->at(i)>0){
        (*M).double_data[0]->at(i) =  (*M).double_data[0]->at(i)/ (*M).double_data[2]->at(i)/pixsize; // rescale mass-weighted values to averages
        (*M).double_data[1]->at(i) =  (*M).double_data[1]->at(i)/ (*M).double_data[2]->at(i)/pixsize;
        (*M).double_data[2]->at(i) =  (*M).double_data[2]->at(i) * MPP; // convert from Npix/steradian to Msun/h / steradian
      }
    }
}

void rescale_map_file_hydro(MapData* M){
    size_t npix = (*M).npixels;
    //TODO: hard-coded for Frontier-E simulation currently
    double MPP = 1342777241.6;
    size_t npix_tot = (*M).npixels_tot;
    double pixsize = (4.*3.141529/npix_tot);
    double xray_scaling = (1./pixsize)/(4.*3.141529);

    // assert that the order is phi, vel, rho
    for (int64_t i=0;i<npix;i++){
      if ((*M).double_data[2]->at(i)>0){
        (*M).double_data[0]->at(i) =  (*M).double_data[0]->at(i)/ (*M).double_data[2]->at(i)/pixsize; // rescale mass-weighted values to averages
        (*M).double_data[1]->at(i) =  (*M).double_data[1]->at(i)/ (*M).double_data[2]->at(i)/pixsize;
        (*M).double_data[7]->at(i) =  (*M).double_data[7]->at(i)/ (*M).double_data[6]->at(i)/xray_scaling; // rescale flux-weighted temperature
      }
    }
}



int redistribute_maps(int numranks, vector<map_properties> maps_send, vector<map_properties>& recv_maps, MapData* M, vector<int> send_count){

    vector<int> send_offset;
    send_offset.resize(numranks);
    send_offset[0] = 0;

    vector<int> recv_offset;
    recv_offset.resize(numranks);
    recv_offset[0] = 0;

    vector<int> recv_count;
    recv_count.resize(numranks);

    // create output for maps 
    MPI_Alltoall( &send_count[0], 1 , MPI_INT, &recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);
    for (int k=1;k<numranks;k++){
        send_offset[k] = send_offset[k-1] + send_count[k-1];
        recv_offset[k] = recv_offset[k-1] + recv_count[k-1];
    }
    recv_maps.resize(recv_offset.back()+recv_count.back());
    //if ((*recv_maps).size()!= recv_offset.back()+recv_count.back()){
    //    if (numranks>1)
    //        (*recv_maps).resize(recv_offset.back()+recv_count.back());
    //}

    // do communication
    if (numranks>1){
        MPI_Alltoallv(&maps_send[0], &send_count[0], &send_offset[0],(*M).map_properties_MPI_Type,\
                  &recv_maps[0], &recv_count[0], &recv_offset[0], (*M).map_properties_MPI_Type, MPI_COMM_WORLD );

        int64_t npix_count = recv_offset.back()+recv_count.back();
        printf("number of pixels on rank = %ld\n",npix_count);

        sort(recv_maps.begin(),recv_maps.end(),comp_pix); // sort by pixel index

        (*M).Resize(npix_count);
        for (int64_t k=0; k<npix_count; k++){
            (*M).Assign(recv_maps.at(k),k);
        }
    }
    else{
	    int64_t npix_count = maps_send.size();
	    sort(maps_send.begin(), maps_send.end(), comp_pix);
	    (*M).Resize(npix_count);
              for (int64_t k=0; k<npix_count; k++){
                (*M).Assign(recv_maps.at(k),k);
            }

    }

    return 0;
}


int compute_ranks( MapData* M, int numranks, vector<int> &send_count, vector<map_properties>& maps_send){
     int64_t npix_per_rank;
     if ((*M).npixels_tot%(int64_t)numranks==0){
        npix_per_rank = (*M).npixels_tot/(int64_t)numranks;
     }
     else{
        npix_per_rank = (*M).npixels_tot/(int64_t)numranks + 1;
     }
     cout << "npix per rank = " << npix_per_rank << endl;
     for (int64_t k=0;k<(*M).npixels;k++){
       map_properties tmp = (*M).GetProperties(k);
       int rank_send = (*M).pix_index[k]/npix_per_rank;
       tmp.rank = rank_send;
       if (rank_send>=numranks){
	       cout << "bad rank send = "<< rank_send << endl;
       }
       assert(rank_send<numranks);
       send_count[rank_send]+=1;
       maps_send.at(k) = tmp;
    }
 return 0;
}



int read_and_redistribute(string file_name, int numranks, MapData* M, vector<map_properties>& maps_send, vector<map_properties>& maps_recv){

    int status; 
    vector<int> send_count(numranks,0);

    read_map_file(M, file_name); // read maps 
    //rescale_map_file(M); //TODO: hard coded for hydro
    rescale_map_file_hydro(M);
    // TURN RESCALING OFF FOR HYDRO
    //
    //
    // if not already the right size let's change it
    //if ((*maps_send).size()!=(*M).npixels){
    maps_send.resize((*M).npixels);
    //}

    status = compute_ranks(M, numranks, send_count, maps_send);
    sort(maps_send.begin(),maps_send.end(),comp_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    status = redistribute_maps(numranks, maps_send, maps_recv, M, send_count);

    return status; 

}

int write_files_slab(MapData &M, int rank, string step_s, vector<float> &out_float){// vector<int> pixel_nums_rank, vector<int64_t> pixel_counts, string zrange_str, string step_s){
        //TODO: currently hard-coded
	// create one file per step 
	int commRank, commRanks;
	MPI_Comm_size(MPI_COMM_WORLD,&commRanks);
        MPI_Comm_rank(MPI_COMM_WORLD,&commRank);

        //string new_name = "/lustre/orion/hep142/proj-shared/prlarsen/maps_hdf5_GO/GO_maps_step_" + step_s + ".hdf5";

	// TODO: uncomment for hydro 
	string new_name = "/lustre/orion/hep142/proj-shared/prlarsen/maps_hdf5_hydro/hydro_maps_step_" + step_s + ".hdf5"; 

	hid_t file_id;
	hid_t file_acctmp; // file access template
	//hid_t dataspace1; // dataspace ids
	//hid_t file_dataspace;
	//hid_t dataset1; 
	hid_t mem_dataspace;
	//hsize_t dims1[SPACE1_RANK] = {SPACE1_DIM1,SPACE1_DIM2}; //TODO fix this
        hsize_t dim_tot[1];
	//dim_tot[0] = 384;
        dim_tot[0] = (hsize_t)(M.npixels_tot); // dimension  TODO: check if this is per rank or total 
        hsize_t dim_chunk[1];
	dim_chunk[0] = 12;

	herr_t ret;
	MPI_Info info = MPI_INFO_NULL;

	file_acctmp = H5Pcreate(H5P_FILE_ACCESS);

        MPI_Barrier(MPI_COMM_WORLD);
        assert(file_acctmp != -1);
	if (commRank==0){
		cout << "HDF5 file creation access succeded" << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
        if (commRank==0){
                cout << "all ranks HDF5 file creation access succeded" << endl;
        }	
        MPI_Barrier(MPI_COMM_WORLD);

	ret = H5Pset_fapl_mpio(file_acctmp, MPI_COMM_WORLD, MPI_INFO_NULL);
        assert(ret != -1);
        if (commRank==0){
                cout << "H5Pset_fapl_mpio succeed" << endl;
        }       

        MPI_Barrier(MPI_COMM_WORLD);


        /* create the file collectively */
        file_id=H5Fcreate(new_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_acctmp);
        assert(file_id != -1);
        if (rank==0){
                cout << "H5Fcreate succeed" << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /* Release file-access template */
        ret=H5Pclose(file_acctmp);
        assert(ret != -1);

        if (rank==0){
                cout << "fat close succeed" << endl;
        }   

	string path_name = "index";
        //dataspace1 = H5Screate_simple (1, dim_tot, NULL);
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        //H5Pset_chunk(plist, 1, dim_chunk);
        //Ha5Pset_deflate(plist, 4);
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank==0){
                cout << "plist created" << endl;
        }
        
	// collective dataset creation
        //dataset1 = H5Dcreate(file_id, "index", H5T_NATIVE_UINT64, dataspace1, H5P_DEFAULT, plist, H5P_DEFAULT);

	MPI_Barrier(MPI_COMM_WORLD);

        //dataset1 = H5Dcreate(file_id, path_name.c_str(), H5T_NATIVE_UINT64, dataspace1, H5P_DEFAULT);
        //assert(dataset1 != -1);
        //if (rank==0){
        //        cout << "dataset1 create succeed" << endl;
        //}   


        uint64_t start_val = M.pix_index[0];
        uint64_t end_val = M.pix_index[M.npixels-1]; // npixels is local 
        //uint64_t start_val_expected = 100663296 * (uint64_t)(commRank);

        hsize_t dim[1];
	//dim[0] = 12;
	dim[0] = (hsize_t)(M.npixels);
	//dim[0] = (int)(M.npixels);

	//assert(M.npixels == 100663296);
        //dim[0] = M.npixels;//end_val - start_val; // local dimension

        hsize_t start[1];
        hsize_t stride[1];
        hsize_t count[1];
        hsize_t block[1];
        //cout << std::numeric_limits<hsize_t>::max() << endl;
	//if (start_val!=100663296*rank){
        //		cout << start_val << " is ne  " << start_val_expected << " on rank " << rank << " =  " <<commRank << endl;//" with max val "<<  std::numeric_limits<hsize_t>::max() <<  endl; 
	//}
	MPI_Barrier(MPI_COMM_WORLD);
	//assert(start_val == start_val_expected);
        start[0] = (hsize_t)start_val;//rank*12;//start_val;
        stride[0] = 1;
	count[0] = (hsize_t)(M.npixels);//12;
        //count[0] = M.npixels;//end_val - start_val; 
        block[0] = 1;

        /*if (rank==0){
                cout << "npixels = " << M.npixels << endl;
		cout << "start val = " << start_val <<endl;
        }*/

        /* create a file dataspace independently */
        //file_dataspace = H5Dget_space (dataset1);
        //assert(file_dataspace != -1);
        //ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
        //assert(ret != -1);
        /* create a memory dataspace independently - 1 dimension, this count*/
        mem_dataspace = H5Screate_simple (1, count, NULL);
        assert (mem_dataspace != -1);
        /* write data independently */
        //ret = H5Dwrite(dataset1, H5T_NATIVE_UINT64, mem_dataspace, file_dataspace, H5P_DEFAULT, &M.pix_index.at(0));

        ///ret = H5Dwrite(dataset1, H5T_NATIVE_UINT64, mem_dataspace, file_dataspace, H5P_DEFAULT, &M.pix_index.at(0));
        //assert(ret != -1);
        //H5Sclose(file_dataspace);

        MPI_Barrier(MPI_COMM_WORLD);
        //if (rank==0){
        //        cout << "first dataset write succeeded  " << endl;
       // }

        out_float.resize(dim[0]);

	for (int i=0; i<N_MAPS;i++){
            hid_t file_dataspace2;
	    hid_t dataspace2;
	    hid_t dataset2;

	    for (size_t j=0; j<dim[0]; j++){
		    out_float[j] = (float)(M.double_data[i]->at(j));
	    }

            path_name = map_names_go[i];
            dataspace2 = H5Screate_simple (1, dim_tot, NULL);
           //dataset1 = H5Dcreate(file_id, "index", H5T_NATIVE_UINT64, dataspace1, H5P_DEFAULT, plist, H5P_DEFAULT);

            dataset2 = H5Dcreate(file_id, path_name.c_str(), H5T_NATIVE_FLOAT, dataspace2, H5P_DEFAULT, plist, H5P_DEFAULT);
            assert(dataset2 != -1);

            /* create a file dataspace independently */
            file_dataspace2 = H5Dget_space (dataset2); // I think this is right
            assert(file_dataspace2 != -1);
            ret=H5Sselect_hyperslab(file_dataspace2, H5S_SELECT_SET, start, stride, count, block);
            assert(ret != -1);
            /* write data independently */
	    // can reuse mem_dataspace
            ret = H5Dwrite(dataset2, H5T_NATIVE_FLOAT, mem_dataspace, file_dataspace2, H5P_DEFAULT, &out_float[0]);
            assert(ret != -1);

	    H5Sclose(file_dataspace2);
	    ret = H5Dclose(dataset2);
	    H5Sclose(dataspace2);

	    assert(ret!=-1);
        }


    H5Pclose(plist);
    H5Sclose(mem_dataspace);
    /* release all IDs created */
    //H5Dclose(dataset1);
    //H5Sclose(dataspace1);
    /* close the file collectively */					    
    H5Fclose(file_id);							    

    // only use rank 0 to write attributes to this file 
     if (rank==0){

    /* uncomment for hydro */

     hid_t atype2, atype3, atype4, atype5, atype6, atype7, atype8, atype9, atype10, atype11;
     hid_t attr_id2, attr_id3, attr_id4, attr_id5, attr_id6, attr_id7, attr_id8, attr_id9, attr_id10, attr_id11;
     hid_t attr2, attr3, attr4, attr5, attr6, attr7, attr8, attr9, attr10, attr11;
     hid_t dataset2, dataset3, dataset4, dataset5, dataset6, dataset7, dataset8, dataset9, dataset10, dataset11;  
					  
     char string2[] = "Average mass-weighted gravitational potential on HEALPix map of nside 16384 in NESTED order";
     char string3[] = "Average mass-weighted line-of-sight velocity in comoving km/s on HEALPix map of nside 16384 in NESTED order";
     char string4[] = "CIC-weighted density map in Msun/h/steradian on HEALPix map of nside 16384 in NESTED order";
     char string5[] = "Total xray flux emission as measured in the ROSAT 0.5keV to 2.0keV band, in units of erg/s/cm^2/steradian on HEALPix map of nside 16384 in NESTED order";
     char string6[] = "Total xray flux emission as measured in the eROSITA soft 0.2keV to 2.3keV band, in units of erg/s/cm^2/steradian on HEALPix map of nside 16384 in NESTED order";
     char string7[] = "Total xray flux emission as measured in the eROSITA hard 2.3keV to 8.0keV band, in units of erg/s/cm^2/steradian on HEALPix map of nside 16384 in NESTED order";
     char string8[] = "Total xray flux emission as measured in a bolometric band covering 0.5keV to 10.0keV, in units of erg/s/cm^2/steradian on HEALPix map of nside 16384 in NESTED order";
     char string9[] = "Average temperature in Kelvin, weighted by the bolometric X-ray flux emission on HEALPix map of nside 16384 in NESTED order";
     char string10[] = "Doppler b-parameter corresponding to the kinetic Sunyaev-Zeldovich effect, unitless, on HEALPix map of nside 16384 in NESTED order";
     char string11[] = "Compton y-parameter corresponding to the thermal Sunyaev-Zeldovich effect, unitless, on HEALPix map of nside 16384 in NESTED order";


     file_id=H5Fopen(new_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
     dataset2=H5Dopen(file_id, "/phi", H5P_DEFAULT);
     dataset3=H5Dopen(file_id, "/vel", H5P_DEFAULT);
     dataset4=H5Dopen(file_id, "/rho", H5P_DEFAULT);
     dataset5=H5Dopen(file_id, "/xray1", H5P_DEFAULT);
     dataset6=H5Dopen(file_id, "/xray2", H5P_DEFAULT);
     dataset7=H5Dopen(file_id, "/xray3", H5P_DEFAULT);
     dataset8=H5Dopen(file_id, "/xray4", H5P_DEFAULT);
     dataset9=H5Dopen(file_id, "/temp", H5P_DEFAULT);
     dataset10=H5Dopen(file_id, "/ksz", H5P_DEFAULT);
     dataset11=H5Dopen(file_id, "/tsz", H5P_DEFAULT);


     attr_id2 = H5Screate(H5S_SCALAR);
     attr_id3 = H5Screate(H5S_SCALAR);
     attr_id4 = H5Screate(H5S_SCALAR);
     attr_id5 = H5Screate(H5S_SCALAR);
     attr_id6 = H5Screate(H5S_SCALAR);
     attr_id7 = H5Screate(H5S_SCALAR);
     attr_id8 = H5Screate(H5S_SCALAR);
     attr_id9 = H5Screate(H5S_SCALAR);
     attr_id10 = H5Screate(H5S_SCALAR);
     attr_id11 = H5Screate(H5S_SCALAR);

     atype2 = H5Tcopy(H5T_C_S1);
     atype3 = H5Tcopy(H5T_C_S1);
     atype4 = H5Tcopy(H5T_C_S1);
     atype5 = H5Tcopy(H5T_C_S1);
     atype6 = H5Tcopy(H5T_C_S1);
     atype7 = H5Tcopy(H5T_C_S1);
     atype8 = H5Tcopy(H5T_C_S1);
     atype9 = H5Tcopy(H5T_C_S1);
     atype10 = H5Tcopy(H5T_C_S1);
     atype11 = H5Tcopy(H5T_C_S1);

     H5Tset_size(atype2, 91);
     H5Tset_size(atype3, 107);
     H5Tset_size(atype4, 90);
     H5Tset_size(atype5, 151);
     H5Tset_size(atype6, 158);
     H5Tset_size(atype7, 158);
     H5Tset_size(atype8, 164);
     H5Tset_size(atype9, 123);
     H5Tset_size(atype10, 130);
     H5Tset_size(atype11, 130);

     attr2 = H5Acreate(dataset2, "Description", atype2, attr_id2, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr2, atype2, string2);
     attr3 = H5Acreate(dataset3, "Description", atype3, attr_id3, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr3, atype3, string3);
     attr4 = H5Acreate(dataset4, "Description", atype4, attr_id4, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr4, atype4, string4);
     attr5 = H5Acreate(dataset5, "Description", atype5, attr_id5, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr5, atype5, string5);
     attr6 = H5Acreate(dataset6, "Description", atype6, attr_id6, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr6, atype6, string6);
     attr7 = H5Acreate(dataset7, "Description", atype7, attr_id7, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr7, atype7, string7);
     attr8 = H5Acreate(dataset8, "Description", atype8, attr_id8, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr8, atype8, string8);
     attr9 = H5Acreate(dataset9, "Description", atype9, attr_id9, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr9, atype9, string9);
     attr10 = H5Acreate(dataset10, "Description", atype10, attr_id10, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr10, atype10, string10);
     attr11 = H5Acreate(dataset11, "Description", atype11, attr_id11, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr11, atype11, string11);

     H5Aclose(attr2);
     H5Aclose(attr3);
     H5Aclose(attr4);
     H5Aclose(attr5);
     H5Aclose(attr6);
     H5Aclose(attr7);
     H5Aclose(attr8);
     H5Aclose(attr9);
     H5Aclose(attr10);
     H5Aclose(attr11);

     H5Sclose(attr_id2);
     H5Sclose(attr_id3);
     H5Sclose(attr_id4);
     H5Sclose(attr_id5);
     H5Sclose(attr_id6);
     H5Sclose(attr_id7);
     H5Sclose(attr_id8);
     H5Sclose(attr_id9);
     H5Sclose(attr_id10);
     H5Sclose(attr_id11);

     H5Dclose(dataset2);
     H5Dclose(dataset3);
     H5Dclose(dataset4);
     H5Dclose(dataset5);
     H5Dclose(dataset6);
     H5Dclose(dataset7);
     H5Dclose(dataset8);
     H5Dclose(dataset9);
     H5Dclose(dataset10);
     H5Dclose(dataset11);

     H5Fclose(file_id);
     /*
 
      // gravity-only, comment out for hydro 
      hid_t atype2, atype3, atype4;
     hid_t attr_id2, attr_id3, attr_id4;
     hid_t attr2, attr3, attr4;
     hid_t dataset2, dataset3, dataset4;  // we need three more datasets defined

     //char string1[] = "HEALPix pixel index corresponding to map of nside 16384 in NESTED order";
     char string2[] = "Average mass-weighted gravitational potential on HEALPix map of nside 16384 in NESTED order";
     char string3[] = "Average mass-weighted line-of-sight velocity in comoving km/s on HEALPix map of nside 16384 in NESTED order";
     char string4[] = "CIC-weighted density map in Msun/h/steradian on HEALPix map of nside 16384 in NESTED order";

     file_id=H5Fopen(new_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
     //dataset1=H5Dopen(file_id, "/index", H5P_DEFAULT);
     dataset2=H5Dopen(file_id, "/phi", H5P_DEFAULT);
     dataset3=H5Dopen(file_id, "/vel", H5P_DEFAULT);
     dataset4=H5Dopen(file_id, "/rho", H5P_DEFAULT);

     //attr_id1 = H5Screate(H5S_SCALAR);
     attr_id2 = H5Screate(H5S_SCALAR);
     attr_id3 = H5Screate(H5S_SCALAR);
     attr_id4 = H5Screate(H5S_SCALAR);

     //atype = H5Tcopy(H5T_C_S1);
     atype2 = H5Tcopy(H5T_C_S1);
     atype3 = H5Tcopy(H5T_C_S1);
     atype4 = H5Tcopy(H5T_C_S1);

     //H5Tset_size(atype, 71);
     H5Tset_size(atype2, 91);
     H5Tset_size(atype3, 107);
     H5Tset_size(atype4, 90);


     //attr1 = H5Acreate(dataset1, "Description", atype, attr_id1, H5P_DEFAULT, H5P_DEFAULT);
     //ret = H5Awrite(attr1, atype, string1);
     attr2 = H5Acreate(dataset2, "Description", atype2, attr_id2, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr2, atype2, string2);
     attr3 = H5Acreate(dataset3, "Description", atype3, attr_id3, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr3, atype3, string3);
     attr4 = H5Acreate(dataset4, "Description", atype4, attr_id4, H5P_DEFAULT, H5P_DEFAULT);
     ret = H5Awrite(attr4, atype4, string4);



     //H5Aclose(attr1);
     H5Aclose(attr2);
     H5Aclose(attr3);
     H5Aclose(attr4);

     //H5Sclose(attr_id1);
     H5Sclose(attr_id2);
     H5Sclose(attr_id3);
     H5Sclose(attr_id4);

     //H5Dclose(dataset1);
     H5Dclose(dataset2);
     H5Dclose(dataset3);
     H5Dclose(dataset4);

     H5Fclose(file_id);

     */
     }

    return 0;
}



