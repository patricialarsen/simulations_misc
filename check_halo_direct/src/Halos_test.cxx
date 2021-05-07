#include <iostream>
#include "PartitionPlus.h"
#include "Halos_test.h"

using namespace std;
using namespace gio;

// Allocate actually allocates memory with optional vector size


void Particles_test::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->num_halos = n;

    this->fof_halo_tag = new vector<int64_t>(this->num_halos);
    this->id = new vector<int64_t>(this->num_halos);

    this->tag2idx = new map<int64_t,int>();

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;
  }
}


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

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;
  }
}

void SODBins_test::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->num_halos = n;

    this->fof_halo_bin_tag = new vector<int64_t>(this->num_halos);
    this->sod_halo_bin = new vector<int32_t>(this->num_halos);


    for (int i=0; i<N_HALO_FLOATS_SOD; ++i)
      this->float_data[i] = new vector<float>(this->num_halos);
    for (int i=0; i<N_HALO_INTS_SOD; ++i)
      this->int_data[i] = new vector<int>(this->num_halos);


    this->tag2idx = new map<int64_t,int>();

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;
  }
}


void Particles_test::Set_MPIType(){
    MPI_Datatype type[3] = { MPI_INT64_T, MPI_INT64_T ,MPI_INT};
    int blocklen[3] = {1,1,1};
    particles_test hp;

    MPI_Aint base;
    MPI_Aint disp[3];

    MPI_Get_address(&hp, &base);
    MPI_Get_address(&hp.fof_halo_tag,     &disp[0]);
    MPI_Get_address(&hp.id,               &disp[1]);
    MPI_Get_address(&hp.rank,             &disp[2]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; 


    MPI_Type_struct(3,blocklen,disp,type,&this->particles_test_MPI_Type);
    MPI_Type_commit(&this->particles_test_MPI_Type);

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


void SODBins_test::Set_MPIType(){
    MPI_Datatype type[5] = { MPI_INT64_T, MPI_INT,  MPI_INT, MPI_FLOAT , MPI_INT};
    int blocklen[5] = {1,1,1,N_HALO_FLOATS_SOD, N_HALO_INTS_SOD};
    sod_binproperties_test hp;

    MPI_Aint base;
    MPI_Aint disp[5];

    MPI_Get_address(&hp, &base);
    MPI_Get_address(&hp.fof_halo_bin_tag,     &disp[0]);
    MPI_Get_address(&hp.sod_halo_bin,   &disp[1]);
    MPI_Get_address(&hp.rank,   &disp[2]);
    MPI_Get_address(&hp.float_data,             &disp[3]);
    MPI_Get_address(&hp.int_data,       &disp[4]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;
    disp[4]-=base;


    MPI_Type_struct(5,blocklen,disp,type,&this->sod_binproperties_MPI_Type);
    MPI_Type_commit(&this->sod_binproperties_MPI_Type);

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

    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}

void SODBins_test::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->num_halos = 0;

    delete this->fof_halo_bin_tag;
    delete this->sod_halo_bin;

    for (int i=0; i<N_HALO_FLOATS_SOD; ++i)
      delete this->float_data[i];
    for (int i=0; i<N_HALO_INTS_SOD; ++i)
      delete this->int_data[i];

    delete this->tag2idx;


    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;
  }
}


void Particles_test::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->num_halos = 0;

    delete this->fof_halo_tag;
    delete this->id;

    delete this->tag2idx;


    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;
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

sod_binproperties_test SODBins_test::GetProperties(size_t idx) {
  if (this->is_allocated) {
    sod_binproperties_test sod_binproperties;

    sod_binproperties.fof_halo_bin_tag =   this->fof_halo_bin_tag->at(idx);
    sod_binproperties.sod_halo_bin = this->sod_halo_bin->at(idx);

    for (int i=0; i<N_HALO_FLOATS_SOD; ++i)
      sod_binproperties.float_data[i] = this->float_data[i]->at(idx);
    for (int i=0; i<N_HALO_INTS_SOD; ++i)
      sod_binproperties.int_data[i] = this->int_data[i]->at(idx);


    return sod_binproperties;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    ;
  }
}


particles_test Particles_test::GetProperties(size_t idx) {
  if (this->is_allocated) {
    particles_test parts;

    parts.fof_halo_tag =   this->fof_halo_tag->at(idx);
    parts.id = this->id->at(idx);


    return parts;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    ;
  }
}



void Halos_test::Erase(size_t i) {
  if (this->is_allocated) {

    this->fof_halo_count->erase(this->fof_halo_count->begin()+i);
    this->fof_halo_tag->erase(this->fof_halo_tag->begin()+i);
    this->sod_halo_count->erase(this->sod_halo_count->begin()+i);
    for (int j=0; j<N_HALO_FLOATS; ++j)
        this->float_data[j]->erase(this->float_data[j]->begin()+i);


    --this->num_halos;

  }
  else {
    cerr << "ERROR: Can't erase -- vectors are not allocated." << endl;
  }
}


void SODBins_test::Erase(size_t i) {
  if (this->is_allocated) {

    this->fof_halo_bin_tag->erase(this->fof_halo_bin_tag->begin()+i);
    this->sod_halo_bin->erase(this->sod_halo_bin->begin()+i);
    for (int j=0; j<N_HALO_FLOATS_SOD; ++j)
        this->float_data[j]->erase(this->float_data[j]->begin()+i);
    for (int j=0; j<N_HALO_INTS_SOD; ++j)
        this->int_data[j]->erase(this->int_data[j]->begin()+i);

    --this->num_halos;

  }
  else {
    cerr << "ERROR: Can't erase -- vectors are not allocated." << endl;
  }
}


void Particles_test::Erase(size_t i) {
  if (this->is_allocated) {

    this->fof_halo_tag->erase(this->fof_halo_tag->begin()+i);
    this->id->erase(this->id->begin()+i);

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
    ;
  }
}

void SODBins_test::PushBack(sod_binproperties_test h) {
  if (this->is_allocated) {
    this->fof_halo_bin_tag->push_back(h.fof_halo_bin_tag);
    this->sod_halo_bin->push_back(h.sod_halo_bin);

    for (int i=0; i<N_HALO_FLOATS_SOD; ++i)
      this->float_data[i]->push_back(h.float_data[i]);
    for (int i=0; i<N_HALO_INTS_SOD; ++i)
      this->int_data[i]->push_back(h.int_data[i]);

    ++this->num_halos;
  }
  else {
    cerr << "ERROR: Can't push back -- vectors are not allocated." << endl;
    ;
  }
}

void Particles_test::PushBack(particles_test h) {
  if (this->is_allocated) {
    this->fof_halo_tag->push_back(h.fof_halo_tag);
    this->id->push_back(h.id);

    ++this->num_halos;
  }
  else {
    cerr << "ERROR: Can't push back -- vectors are not allocated." << endl;
    ;
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
    ;
  }
}
void SODBins_test::Resize(size_t n) {
  if (this->is_allocated) {
    this->num_halos = n;
    this->fof_halo_bin_tag->resize(this->num_halos);
    this->sod_halo_bin->resize(this->num_halos);

    for (int i=0; i<N_HALO_FLOATS_SOD; ++i)
      this->float_data[i]->resize(this->num_halos);
    for (int i=0; i<N_HALO_INTS_SOD; ++i)
      this->int_data[i]->resize(this->num_halos);
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}

void Particles_test::Resize(size_t n) {
  if (this->is_allocated) {
    this->num_halos = n;
    this->fof_halo_tag->resize(this->num_halos);
    this->id->resize(this->num_halos);

  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}

