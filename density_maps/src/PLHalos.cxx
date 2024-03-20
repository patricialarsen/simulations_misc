#include <iostream>
#include "PLHalos.h"

using namespace std;

void PLHalos::Allocate(size_t n) {
  if (!this->is_allocated) {

    this->is_allocated = true;
    this->nparticles = n;

    this->pix_index = new vector<int>(this->nparticles);
    this-> rank = new vector<int>(this->nparticles);

    
    //cout << "N_FLOATS "<< N_FLOATS;
    for (int i=0; i<N_FLOATS_HALOS; ++i)
      this->float_data[i] = new vector<float>(this->nparticles);

    for (int i=0; i<N_DOUBLES_HALOS; ++i)
      this->double_data[i] = new vector<double>(this->nparticles);

    for (int i=0; i<N_INTS_HALOS; ++i)
      this->int_data[i] = new vector<int>(this->nparticles);

    for (int i=0; i<N_INT64S_HALOS; ++i)
      this->int64_data[i] = new vector<int64_t>(this->nparticles);

    for (int i=0; i<N_MASKS_HALOS; ++i)
      this->mask_data[i] = new vector<uint16_t>(this->nparticles);

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;
  }
}



void PLHalos::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->nparticles = 0;

    delete this->pix_index;
    delete this-> rank;

    for (int i=0; i<N_FLOATS_HALOS; ++i)
      delete this->float_data[i];
    for (int i=0; i<N_DOUBLES_HALOS; ++i)
      delete this->double_data[i];
    for (int i=0; i<N_INTS_HALOS; ++i)
      delete this->int_data[i];
    for (int i=0; i<N_INT64S_HALOS; ++i)
      delete this->int64_data[i];
    for (int i=0; i<N_MASKS_HALOS; ++i)
      delete this->mask_data[i];
    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
  }
}


void PLHalos::Resize(size_t n) {
  if (this->is_allocated) {
    this->nparticles = n;

    this->pix_index->resize(this->nparticles);
    this->rank->resize(this->nparticles);

    for (int i=0; i<N_FLOATS_HALOS; ++i)
      this->float_data[i]->resize(this->nparticles);
    for (int i=0; i<N_DOUBLES_HALOS; ++i)
      this->double_data[i]->resize(this->nparticles);
    for (int i=0; i<N_INTS_HALOS; ++i)
      this->int_data[i]->resize(this->nparticles);
    for (int i=0; i<N_INT64S_HALOS; ++i)
      this->int64_data[i]->resize(this->nparticles);
    for (int i=0; i<N_MASKS_HALOS; ++i)
      this->mask_data[i]->resize(this->nparticles);
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}

