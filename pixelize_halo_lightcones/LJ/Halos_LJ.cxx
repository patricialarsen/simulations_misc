#include <iostream>
//#include "PartitionPlus.h"
#include "Halos_LJ.h"

using namespace std;
using namespace gio;

// Allocate actually allocates memory with optional vector size

void Halos_test::Allocate(size_t n) {
  if (!this->is_allocated) {
    this->is_allocated = true;
    this->num_halos = n;

    this->pix_index = new vector<int>(this->num_halos);
    this->fof_halo_tag = new vector<int64_t>(this->num_halos);
    this-> rank = new vector<int>(this->num_halos);

    for (int i=0; i<N_HALO_FLOATS; ++i)
      this->float_data[i] = new vector<float>(this->num_halos);

    for (int i=0; i<N_LIGHTCONE_FLOATS; ++i)
      this->lightcone_data[i] = new vector<float>(this->num_halos);

    for (int i=0; i<N_HALO_FLOATS_E; ++i)
      this->ellipticity_data[i] = new vector<float>(this->num_halos);

    //this->tag2idx = new map<int64_t,int>(); 

  }
  else {
    cerr << "ERROR: Can't allocate new vectors --\
        vectors are already allocated." << endl;
    ;
  }
}


void Halos_test::Set_MPIType(){
    MPI_Datatype type[6] = { MPI_INT64_T, MPI_INT, MPI_INT, MPI_FLOAT ,MPI_FLOAT, MPI_FLOAT};
    int blocklen[6] = {1, 1, 1, N_HALO_FLOATS, N_LIGHTCONE_FLOATS, N_HALO_FLOATS_E};
    halo_properties hp;

    MPI_Aint base;
    MPI_Aint disp[6];

    MPI_Get_address(&hp, &base);
    MPI_Get_address(&hp.fof_halo_tag,     &disp[0]);
    MPI_Get_address(&hp.rank,             &disp[1]);
    MPI_Get_address(&hp.pix_index,        &disp[2]);
    MPI_Get_address(&hp.float_data,       &disp[3]);
    MPI_Get_address(&hp.lightcone_data,       &disp[4]);
    MPI_Get_address(&hp.ellipticity_data, &disp[5]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;
    disp[4]-=base; disp[5]-=base;

    MPI_Type_struct(6,blocklen,disp,type,&this->halo_properties_MPI_Type);
    MPI_Type_commit(&this->halo_properties_MPI_Type);

}



void Halos_test::Deallocate() {

  if (this->is_allocated) {
    this->is_allocated = false;
    this->num_halos = 0;

    delete this->pix_index;
    delete this->fof_halo_tag;
    delete this-> rank;

    for (int i=0; i<N_HALO_FLOATS; ++i)
      delete this->float_data[i];
    for (int i=0; i<N_LIGHTCONE_FLOATS; ++i)
      delete this->lightcone_data[i];
    for (int i=0; i<N_HALO_FLOATS_E; ++i)
      delete this->ellipticity_data[i];


    //delete this->tag2idx;

    this->is_allocated=false;
  }
  else {
    cerr << "ERROR: Can't deallocate -- vectors are not allocated."\
        << endl;
    ;// do some graceful exit of the program
  }
}




halo_properties Halos_test::GetProperties(size_t idx) {
  if (this->is_allocated) {
    halo_properties halo_properties;

    halo_properties.fof_halo_tag =   this->fof_halo_tag->at(idx);

    halo_properties.rank = this->rank->at(idx);

    halo_properties.pix_index = this->pix_index->at(idx);

    for (int i=0; i<N_HALO_FLOATS; ++i)
      halo_properties.float_data[i] = this->float_data[i]->at(idx);

    for (int i=0; i<N_LIGHTCONE_FLOATS; ++i)
      halo_properties.lightcone_data[i] = this->lightcone_data[i]->at(idx);
 
   
    for (int i=0; i<N_HALO_FLOATS_E; ++i)
      halo_properties.ellipticity_data[i] = this->ellipticity_data[i]->at(idx);



    return halo_properties;
  }
  else {
    cerr << "ERROR: Can't get properties -- vectors are not allocated."\
        << endl;
    halo_properties halo_properties;
    return halo_properties;
    ;// do some graceful exit of the program
  }
}



void Halos_test::Erase(size_t i) {
  if (this->is_allocated) {

    this->rank->erase(this->rank->begin()+i);
    this->pix_index->erase(this->pix_index->begin()+i);
    this->fof_halo_tag->erase(this->fof_halo_tag->begin()+i);
    for (int j=0; j<N_HALO_FLOATS; ++j)
        this->float_data[j]->erase(this->float_data[j]->begin()+i);
    for (int j=0; j<N_LIGHTCONE_FLOATS; ++j)
        this->lightcone_data[j]->erase(this->lightcone_data[j]->begin()+i);
    for (int j=0; j<N_HALO_FLOATS_E; ++j)
        this->ellipticity_data[j]->erase(this->ellipticity_data[j]->begin()+i);

    --this->num_halos;

  }
  else {
    cerr << "ERROR: Can't erase -- vectors are not allocated." << endl;
  }
}



void Halos_test::PushBack(halo_properties h) {
  if (this->is_allocated) {
    this->pix_index->push_back(h.pix_index);
    this->fof_halo_tag->push_back(h.fof_halo_tag);
    this->rank->push_back(h.rank);
    for (int i=0; i<N_HALO_FLOATS; ++i)
      this->float_data[i]->push_back(h.float_data[i]);
    for (int i=0; i<N_LIGHTCONE_FLOATS; ++i)
      this->lightcone_data[i]->push_back(h.lightcone_data[i]);
    for (int i=0; i<N_HALO_FLOATS_E; ++i)
      this->ellipticity_data[i]->push_back(h.ellipticity_data[i]);

    ++this->num_halos;
  }
  else {
    cerr << "ERROR: Can't push back -- vectors are not allocated." << endl;
    ;
  }
}



void Halos_test::Assign(halo_properties h, size_t j){
    if (this->is_allocated){
      this->pix_index->at(j) = h.pix_index;
      this->fof_halo_tag->at(j) = h.fof_halo_tag;
      this->rank->at(j) = h.rank;
    for (int i=0; i<N_HALO_FLOATS; ++i)
      this->float_data[i]->at(j)=h.float_data[i];
    for (int i=0; i<N_LIGHTCONE_FLOATS; ++i)
      this->lightcone_data[i]->at(j)=h.lightcone_data[i];
    for (int i=0; i<N_HALO_FLOATS_E; ++i)
      this->ellipticity_data[i]->at(j)=h.ellipticity_data[i];
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}




// use Resize(0) as a substitute for std::vector clear()
void Halos_test::Resize(size_t n) {
  if (this->is_allocated) {
    this->num_halos = n;
    this->pix_index->resize(this->num_halos);
    this->fof_halo_tag->resize(this->num_halos);
    this->rank->resize(this->num_halos);
    for (int i=0; i<N_HALO_FLOATS; ++i)
      this->float_data[i]->resize(this->num_halos);

    for (int i=0; i<N_LIGHTCONE_FLOATS; ++i)
      this->lightcone_data[i]->resize(this->num_halos);
  
    for (int i=0; i<N_HALO_FLOATS_E;++i)
      this->ellipticity_data[i]->resize(this->num_halos);
  }
  else {
    cerr << "ERROR: Can't resize -- vectors are not allocated." << endl;
    ;
  }
}

