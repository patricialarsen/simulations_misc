#include <iostream>
#include "PartitionPlus.h"
#include "SZ_class.h"

using namespace std;
using namespace gio;

// Allocate actually allocates memory with optional vector size
void SZ_class::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->num_parts = n;

    this->pix_num = new vector<int64_t>(this->num_parts);
    this->rank = new vector<int>(this->num_parts); 
    this->tsz = new vector<double>(this->num_parts);
    this->ksz = new vector<double>(this->num_parts);
  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;// do some graceful exit of the program
  }
}

void SZ_class::Set_MPIType(){
    MPI_Datatype type[4] = { MPI_INT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };
    int blocklen[4] = {1,1,1,1};

    sz_props sz;

    MPI_Aint base;
    MPI_Aint disp[4];

    MPI_Get_address(&sz, &base);
    MPI_Get_address(&sz.pix_num,  &disp[0]);
    MPI_Get_address(&sz.tsz,   &disp[1]);
    MPI_Get_address(&sz.ksz,    &disp[2]);
    MPI_Get_address(&sz.rank,    &disp[3]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;

    MPI_Type_struct(4,blocklen,disp,type,&this->sz_properties_MPI_Type);
    MPI_Type_commit(&this->sz_properties_MPI_Type);

}

void SZ_class::Deallocate() {

  if (this->is_allocated) {

    this->num_parts = 0;

    delete this-> pix_num;
    delete this-> rank;
    delete this-> tsz;
    delete this-> ksz;
    
    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}

sz_props SZ_class::GetProperties(size_t idx) {
  if (this->is_allocated) {
    sz_props sz_properties;

    sz_properties.rank =   this->rank->at(idx);
    sz_properties.pix_num =   this->pix_num->at(idx);
    sz_properties.tsz =   this->tsz->at(idx);
    sz_properties.ksz =   this->ksz->at(idx);

    return sz_properties;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}


// use Resize(0) as a substitute for std::vector clear()
void SZ_class::Resize(size_t n) {
  if (this->is_allocated) {
    this->num_parts = n;
    this->pix_num->resize(this->num_parts);
    this->rank->resize(this->num_parts);
    this->tsz->resize(this->num_parts);
    this->ksz->resize(this->num_parts);
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;// do some graceful exit of the program
  }
}

