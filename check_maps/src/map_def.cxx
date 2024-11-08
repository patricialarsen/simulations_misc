#include <iostream>
#include "map_def.h"

using namespace std;


void MapData::Resize(size_t n) {
  if (this->is_allocated) {
    this->npixels = n;
    this->pix_index.resize(this->npixels);
    for (int i=0; i<N_MAPS; ++i)
      this->double_data[i]->resize(this->npixels);
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
      this->double_data[i]->resize(n);
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
      this->double_data[i] = new vector<double>(this->npixels);

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
      this->double_data[i] = new vector<double>(this->npixels);

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

