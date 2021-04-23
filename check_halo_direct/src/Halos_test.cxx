#include <iostream>
#include "PartitionPlus.h"
#include "Halos_test.h"

using namespace std;
using namespace gio;

// Allocate actually allocates memory with optional vector size
void Halos_test::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->num_halos = n;

    this->fof_halo_tag = new vector<int64_t>(this->num_halos);
    this->fof_halo_count = new vector<int32_t>(this->num_halos); 
    this->sod_halo_count = new vector<int64_t>(this->num_halos);


    for (int i=0; i<N_HALO_FLOATS; ++i)
      this->float_data[i] = new vector<float>(this->num_halos);

    this->tag2idx = new map<int64_t,int>(); 

    this->pid = new allocated_vector<int64_t>(); 
    this->tag = new allocated_vector<int64_t>(); 
  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;// do some graceful exit of the program
  }
}

void Halos_test::Set_MPIType(){
    MPI_Datatype type[5] = { MPI_INT64_T, MPI_INT, MPI_INT64_T, MPI_INT, MPI_FLOAT };
    int blocklen[5] = {1,1,1,1,N_HALO_FLOATS};
    halo_properties_test hp;

    MPI_Aint base;
    MPI_Aint disp[5];

    MPI_Get_address(&hp, &base);
    MPI_Get_address(&hp.fof_halo_tag,     &disp[0]);
    MPI_Get_address(&hp.fof_halo_count,   &disp[1]);
    MPI_Get_address(&hp.sod_halo_count,   &disp[2]);
    MPI_Get_address(&hp.rank,             &disp[3]);
    MPI_Get_address(&hp.float_data,       &disp[4]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;
    disp[4]-=base;


    MPI_Type_struct(5,blocklen,disp,type,&this->halo_properties_MPI_Type);
    MPI_Type_commit(&this->halo_properties_MPI_Type);

}


void Halos_test::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->num_halos = 0;

    delete this->fof_halo_tag;
    delete this->fof_halo_count;
    delete this->sod_halo_count;

    for (int i=0; i<N_HALO_FLOATS; ++i)
      delete this->float_data[i];

    delete this->tag2idx;

    delete this->pid;
    delete this->tag;

    if (this->buffer!=NULL)
      free(this->buffer);

    this->buffer=NULL;
    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}

halo_properties_test Halos_test::GetProperties(size_t idx) {
  if (this->is_allocated) {
    halo_properties_test halo_properties;

    halo_properties.fof_halo_tag =   this->fof_halo_tag->at(idx);
    halo_properties.fof_halo_count = this->fof_halo_count->at(idx);
    halo_properties.sod_halo_count = this->sod_halo_count->at(idx);


    for (int i=0; i<N_HALO_FLOATS; ++i)
      halo_properties.float_data[i] = this->float_data[i]->at(idx);
    
    return halo_properties;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}

void Halos_test::Erase(size_t i) {
  if (this->is_allocated) {

    this->fof_halo_count->erase(this->fof_halo_count->begin()+i);
    this->fof_halo_tag->erase(this->fof_halo_tag->begin()+i);
    this->sod_halo_count->erase(this->sod_halo_count->begin()+i);
    for (int i=0; i<N_HALO_FLOATS; ++i)
        this->float_data[i]->erase(this->float_data[i]->begin()+i);


    --this->num_halos;

  }
  else {
    cerr << "ERROR: Can't erase -- vectors are not allocated." << endl;
  }
}



void Halos_test::PushBack(halo_properties_test h) {
  if (this->is_allocated) {
    this->fof_halo_tag->push_back(h.fof_halo_tag);
    this->fof_halo_count->push_back(h.fof_halo_count);
    this->sod_halo_count->push_back(h.sod_halo_count);

    for (int i=0; i<N_HALO_FLOATS; ++i)
      this->float_data[i]->push_back(h.float_data[i]);

    ++this->num_halos;
  }
  else {
    cerr << "ERROR: Can't push back -- vectors are not allocated." << endl;
    ;// do some graceful exit of the program
  }
}

// use Resize(0) as a substitute for std::vector clear()
void Halos_test::Resize(size_t n) {
  if (this->is_allocated) {
    this->num_halos = n;
    this->fof_halo_count->resize(this->num_halos);
    this->fof_halo_tag->resize(this->num_halos);
    this->sod_halo_count->resize(this->num_halos);

    for (int i=0; i<N_HALO_FLOATS; ++i)
      this->float_data[i]->resize(this->num_halos);
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;// do some graceful exit of the program
  }
}

