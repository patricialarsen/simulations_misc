#include <iostream>
#include "PartitionPlus.h"
#include "LC_test.h"

using namespace std;
using namespace gio;

// Allocate actually allocates memory with optional vector size
void LC_test::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->num_parts = n;

    this->id = new vector<int64_t>(this->num_parts);
    this->replication = new vector<int32_t>(this->num_parts); 


    for (int i=0; i<N_LC_FLOATS; ++i)
      this->float_data[i] = new vector<float>(this->num_parts);

    this->tag2idx = new map<int64_t,int>(); 

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;// do some graceful exit of the program
  }
}

void LC_test::Set_MPIType(){
    MPI_Datatype type[4] = { MPI_INT64_T, MPI_INT, MPI_INT, MPI_FLOAT };
    int blocklen[4] = {1,1,1,N_LC_FLOATS};
    lc_properties_test hp;

    MPI_Aint base;
    MPI_Aint disp[4];

    MPI_Get_address(&hp, &base);
    MPI_Get_address(&hp.id,     &disp[0]);
    MPI_Get_address(&hp.replication,   &disp[1]);
    MPI_Get_address(&hp.rank,             &disp[2]);
    MPI_Get_address(&hp.float_data,       &disp[3]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;


    MPI_Type_struct(4,blocklen,disp,type,&this->LC_properties_MPI_Type);
    MPI_Type_commit(&this->LC_properties_MPI_Type);

}


void LC_test::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->num_parts = 0;

    delete this->id;
    delete this->replication;

    for (int i=0; i<N_LC_FLOATS; ++i)
      delete this->float_data[i];

    delete this->tag2idx;

    //if (this->buffer!=NULL)
    //  free(this->buffer);

   // this->buffer=NULL;
    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}

lc_properties_test LC_test::GetProperties(size_t idx) {
  if (this->is_allocated) {
    lc_properties_test lc_properties;

    lc_properties.id =   this->id->at(idx);
    lc_properties.replication = this->replication->at(idx);

    for (int i=0; i<N_LC_FLOATS; ++i)
      lc_properties.float_data[i] = this->float_data[i]->at(idx);
    
    return lc_properties;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}

void LC_test::Erase(size_t i) {
  if (this->is_allocated) {

    this->id->erase(this->id->begin()+i);
    this->replication->erase(this->replication->begin()+i);
    for (int i=0; i<N_LC_FLOATS; ++i)
        this->float_data[i]->erase(this->float_data[i]->begin()+i);

    --this->num_parts;

  }
  else {
    cerr << "ERROR: Can't erase -- vectors are not allocated." << endl;
  }
}



void LC_test::PushBack(lc_properties_test h) {
  if (this->is_allocated) {
    this->id->push_back(h.id);
    this->replication->push_back(h.replication);

    for (int i=0; i<N_LC_FLOATS; ++i)
      this->float_data[i]->push_back(h.float_data[i]);

    ++this->num_parts;
  }
  else {
    cerr << "ERROR: Can't push back -- vectors are not allocated." << endl;
    ;// do some graceful exit of the program
  }
}

// use Resize(0) as a substitute for std::vector clear()
void LC_test::Resize(size_t n) {
  if (this->is_allocated) {
    this->num_parts = n;
    this->id->resize(this->num_parts);
    this->replication->resize(this->num_parts);

    for (int i=0; i<N_LC_FLOATS; ++i)
      this->float_data[i]->resize(this->num_parts);
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;// do some graceful exit of the program
  }
}

