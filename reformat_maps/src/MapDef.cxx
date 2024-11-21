#include <iostream>
//#include "map_def.h"
#include "MapDef.h"

using namespace std;

void MapData::Resize(size_t n) {
  if (this->is_allocated) {
    this->npixels = n;
    this->pix_index.resize(n);
    for (int i=0; i<N_MAPS; ++i)
      (*(this->double_data[i])).resize(n); // check this
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;// do some graceful exit of the program
  }
}

void MapDataHydro::Resize(size_t n) {
  if (this->is_allocated) {
    this->npixels = n;
    this->pix_index.resize(n);
    for (int i=0; i<N_MAPS_HYDRO; ++i)
      (*(this->double_data[i])).resize(n); // check this
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;// do some graceful exit of the program
  }
}


void MapData::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->npixels = 0;
    this->npixels_tot = 0;
    this->pix_index.resize(0);
    for (int i=0; i<N_MAPS; ++i)
      delete this->double_data[i];
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}

void MapDataHydro::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->npixels = 0;
    this->npixels_tot = 0;
    this->pix_index.resize(0);
    for (int i=0; i<N_MAPS_HYDRO; ++i)
      delete this->double_data[i];
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}


void MapData::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->npixels = n;
    this->pix_index.resize(n);
    for (int i=0; i<N_MAPS; ++i)
      this->double_data[i] = new vector<double>(n);

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;// do some graceful exit of the program
  }
}

void MapDataHydro::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->npixels = n;
    this->pix_index.resize(n);

    for (int i=0; i<N_MAPS_HYDRO; ++i)
      this->double_data[i] = new vector<double>(n);

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;// do some graceful exit of the program
  }
}

void MapData::Clear(size_t n) {
    for (int i=0; i<N_MAPS; ++i){
      for (int64_t jj=0; jj<n; jj++){
        this->double_data[i]->at(jj) = 0;
      }
    }
}


void MapDataHydro::Clear(size_t n) {
    for (int i=0; i<N_MAPS_HYDRO; ++i){
      for (int64_t jj=0; jj<n; jj++){
        this->double_data[i]->at(jj) = 0;
      }
    }
}
//
// Edit the below depending on what is needed
//
//
void MapData::Set_MPIType(){
	MPI_Datatype type[3] = {MPI_INT64_T, MPI_INT, MPI_DOUBLE};
	int blocklen[3] = {1,1,N_MAPS};
	map_properties mp;

	MPI_Aint base;
	MPI_Aint disp[3];
        MPI_Get_address(&mp, &base);
        MPI_Get_address(&mp.rank,             &disp[1]);
        MPI_Get_address(&mp.pix_index,        &disp[0]);
        MPI_Get_address(&mp.map_data,      &disp[2]);

        disp[0]-=base; disp[1]-=base; disp[2]-=base;

    MPI_Type_struct(3,blocklen,disp,type,&this->map_properties_MPI_Type);
    MPI_Type_commit(&this->map_properties_MPI_Type);

}
void MapDataHydro::Set_MPIType(){
        MPI_Datatype type[3] = {MPI_INT64_T, MPI_INT, MPI_DOUBLE};
        int blocklen[3] = {1,1,N_MAPS_HYDRO};
        map_properties_hydro mp;

        MPI_Aint base;
        MPI_Aint disp[3];
        MPI_Get_address(&mp, &base);
        MPI_Get_address(&mp.rank,             &disp[1]);
        MPI_Get_address(&mp.pix_index,        &disp[0]);
        MPI_Get_address(&mp.map_data,      &disp[2]);

        disp[0]-=base; disp[1]-=base; disp[2]-=base;

    MPI_Type_struct(3,blocklen,disp,type,&this->map_properties_hydro_MPI_Type);
    MPI_Type_commit(&this->map_properties_hydro_MPI_Type);

}

map_properties MapData::GetProperties(size_t i) {
  if (this->is_allocated) {

    map_properties map_properties;

    map_properties.pix_index = this->pix_index.at(i); // pix_index is a vector
    
    for (int j=0; j<N_MAPS; j++)
	    map_properties.map_data[j] = this->double_data[j]->at(i);

    return map_properties;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    map_properties map_properties;
    return map_properties;
    ;// do some graceful exit of the program
  }
}

map_properties_hydro MapDataHydro::GetProperties(size_t i) {
  if (this->is_allocated) {

    map_properties_hydro map_properties;

    map_properties.pix_index = this->pix_index.at(i);

    for (int j=0; j<N_MAPS_HYDRO; j++)
            map_properties.map_data[j] = this->double_data[j]->at(i);

    return map_properties;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    map_properties_hydro map_properties;
    return map_properties;
    ;// do some graceful exit of the program
  }
}


void MapData::Assign(map_properties m, size_t j){
    if (this->is_allocated){
      this->pix_index.at(j) = m.pix_index;
    for (int i=0; i<N_MAPS; ++i)
      this->double_data[i]->at(j)=m.map_data[i];
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}


void MapDataHydro::Assign(map_properties_hydro m, size_t j){
    if (this->is_allocated){
      this->pix_index.at(j) = m.pix_index;
    for (int i=0; i<N_MAPS_HYDRO; ++i)
      this->double_data[i]->at(j)=m.map_data[i];
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}


