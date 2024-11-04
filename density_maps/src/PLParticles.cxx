#include <iostream>
#include "PLParticles.h"

using namespace std;

void PLParticles::Allocate(size_t n) {
  if (!this->is_allocated) {

    this->is_allocated = true;
    this->nparticles = n;

    this->pix_index = new vector<int>(this->nparticles);
    this-> rank = new vector<int>(this->nparticles);

    
    //cout << "N_FLOATS "<< N_FLOATS;
    for (int i=0; i<N_FLOATS; ++i)
      this->float_data[i] = new vector<float>(this->nparticles);

    for (int i=0; i<N_DOUBLES; ++i)
      this->double_data[i] = new vector<double>(this->nparticles);

    for (int i=0; i<N_INTS; ++i)
      this->int_data[i] = new vector<int>(this->nparticles);

    for (int i=0; i<N_INT64S; ++i)
      this->int64_data[i] = new vector<int64_t>(this->nparticles);

    for (int i=0; i<N_MASKS; ++i)
      this->mask_data[i] = new vector<uint16_t>(this->nparticles);

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;
  }
}



void PLParticles::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->nparticles = 0;

    delete this->pix_index;
    delete this-> rank;

    for (int i=0; i<N_FLOATS; ++i)
      delete this->float_data[i];
    for (int i=0; i<N_DOUBLES; ++i)
      delete this->double_data[i];
    for (int i=0; i<N_INTS; ++i)
      delete this->int_data[i];
    for (int i=0; i<N_INT64S; ++i)
      delete this->int64_data[i];
    for (int i=0; i<N_MASKS; ++i)
      delete this->mask_data[i];
    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
  }
}


void PLParticles::Resize(size_t n) {
  if (this->is_allocated) {
    this->nparticles = n;

    this->pix_index->resize(this->nparticles);
    this->rank->resize(this->nparticles);

    for (int i=0; i<N_FLOATS; ++i)
      this->float_data[i]->resize(this->nparticles);
    for (int i=0; i<N_DOUBLES; ++i)
      this->double_data[i]->resize(this->nparticles);
    for (int i=0; i<N_INTS; ++i)
      this->int_data[i]->resize(this->nparticles);
    for (int i=0; i<N_INT64S; ++i)
      this->int64_data[i]->resize(this->nparticles);
    for (int i=0; i<N_MASKS; ++i)
      this->mask_data[i]->resize(this->nparticles);
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}

void PLParticles::Transform(const std::vector<int>& indices) {
    if (this->is_allocated) {
        this->nparticles = indices.size();

        for (size_t i = 0; i < indices.size(); ++i) {//I dont think this is thread safe so not using pragma
            this->pix_index->at(i) = this->pix_index->at(indices[i]);
            this->rank->at(i) = this->rank->at(indices[i]);
            for (int j = 0; j < N_FLOATS; ++j)
                this->float_data[j]->at(i) = this->float_data[j]->at(indices[i]);
            for (int j = 0; j < N_DOUBLES; ++j)
                this->double_data[j]->at(i) = this->double_data[j]->at(indices[i]);
            for (int j = 0; j < N_INTS; ++j)
                this->int_data[j]->at(i) = this->int_data[j]->at(indices[i]);
            for (int j = 0; j < N_INT64S; ++j)
                this->int64_data[j]->at(i) = this->int64_data[j]->at(indices[i]);
            for (int j = 0; j < N_MASKS; ++j)
                this->mask_data[j]->at(i) = this->mask_data[j]->at(indices[i]);
        }

        // Resize the data
        this->pix_index->resize(this->nparticles);
        this->rank->resize(this->nparticles);

        for (int i = 0; i < N_FLOATS; ++i)
            this->float_data[i]->resize(this->nparticles);
        for (int i = 0; i < N_DOUBLES; ++i)
            this->double_data[i]->resize(this->nparticles);
        for (int i = 0; i < N_INTS; ++i)
            this->int_data[i]->resize(this->nparticles);
        for (int i = 0; i < N_INT64S; ++i)
            this->int64_data[i]->resize(this->nparticles);
        for (int i = 0; i < N_MASKS; ++i)
            this->mask_data[i]->resize(this->nparticles);
    } else {
        std::cerr << "ERROR: Can't resize -- vectors are not allocated." << std::endl;
    }
}

