#include <iostream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <stdlib.h>

#include "GetRusageMPI.h"
#include "PartitionPlus.h"
#include "HaloDistribute_mt.h"
#include "SimpleTimings.h"

#define GB_BYTES 1e9
#define MB_BYTES 1e6

#define MSG_SIZE 100*MB_BYTES

using namespace std;
using namespace gio;

inline int64_t gen_tree_node_index(int idx, int rank, int timestep) {

  int64_t nodeindex = 0;

  if (idx==-1)
    return -1;

  assert(idx >= 0);
  assert(rank >= 0 && rank < cosmotk::PartitionPlus::getNumProc());
  assert(timestep >= 0);

  nodeindex |= ((int64_t)timestep << 52);
  nodeindex |= ((int64_t)rank     << 32);
  nodeindex |= ((int64_t)idx);
  return nodeindex;
}

inline bool dis_tree_node_index(const int64_t tni, int &idx, int &rank, int &timestep) {
  if (tni==-1) {
    idx = -1; rank = cosmotk::PartitionPlus::getMyProc(); timestep = -1;
    return false;
  }

  idx =  (0x00000000FFFFFFFF & tni);        // idx is lower 32 bits
  rank = (0x000FFFFF00000000 & tni) >> 32;  // 20 bit rank
  timestep = tni >> 52;                     // 12 bit timestep

  return true;
}

bool sort_by_pid(const pid_tag_t &a, const pid_tag_t &b) {
  return a.pid < b.pid;
}

bool sort_by_dest(const pid_tag_t &a, const pid_tag_t &b) {
  return a.rank< b.rank;
}

bool sort_halos_by_dest(const halo_properties_t &a, const halo_properties_t &b) {
  return a.rank < b.rank;
}

inline void periodic_boundaries(vector<POSVEL_T>* x, vector<POSVEL_T>* y, vector<POSVEL_T>* z,\
    POSVEL_T boxSize) {

  assert(x->size()==y->size());
  assert(x->size()==z->size());

  size_t count = x->size();

  for (size_t i = 0; i < count; ++i) {
    if (x->at(i) >= boxSize)
      x->at(i) -= boxSize;
    if (y->at(i) >= boxSize)
      y->at(i) -= boxSize;
    if (z->at(i) >= boxSize)
      z->at(i) -= boxSize;

    if (x->at(i) < 0.0)
      x->at(i) += boxSize;
    if (y->at(i) < 0.0)
      y->at(i) += boxSize;
    if (z->at(i) < 0.0)
      z->at(i) += boxSize;
  }
}

namespace cosmotk {

// constructor needs pointers to halo memory
HaloDistribute::HaloDistribute(Halos* halos) : halos(halos) {
  ;
}

// destructor
HaloDistribute::~HaloDistribute() {
  ;
}

// Set parameters for halo distribution
void HaloDistribute::set_parameters(
    bool redistribute,
    float pad_factor,
    POSVEL_T rL,
    float halo_loss_threshold) {

  this->redistribute = redistribute;
  this->pad_factor = pad_factor;
  this->boxSize = rL;
  this->halo_loss_threshold = halo_loss_threshold;
 
  if (redistribute)
    this->MB = GenericIO::MismatchRedistribute;
  else
    this->MB = GenericIO::MismatchDisallowed;
}

void HaloDistribute::initialize() {

  this->numProc = PartitionPlus::getNumProc();
  this->myProc =  PartitionPlus::getMyProc();

  // Get the number of ranks in each dimension
  PartitionPlus::getDecompSize(this->layoutSize);

  // Get my position within the Cartesian topology
  PartitionPlus::getMyPosition(this->layoutPos);

  // Get neighbor of this rank
  PartitionPlus::getNeighbors(this->neighbor);

  // MPI Type for halo_properties_t
  {
    halo_properties_t hp;
    MPI_Aint base;
    MPI_Aint disp[9];

    MPI_Get_address(&hp, &base);
    MPI_Get_address(&hp.fof_halo_tag,     &disp[0]);
    MPI_Get_address(&hp.fof_halo_count,   &disp[1]);
    MPI_Get_address(&hp.sod_halo_count,   &disp[2]);
    MPI_Get_address(&hp.rank,             &disp[3]);
    MPI_Get_address(&hp.descendant,       &disp[4]);
    MPI_Get_address(&hp.tree_node_mass,   &disp[5]);
    MPI_Get_address(&hp.float_data,       &disp[6]);
    MPI_Get_address(&hp.tree_node_index,  &disp[7]);
    MPI_Get_address(&hp.desc_node_index,  &disp[8]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base; disp[3]-=base;
    disp[4]-=base; disp[5]-=base; disp[6]-=base; disp[7]-=base;
    disp[8]-=base;

    int blocklen[9] = {1, 1, 1, 1, 1, 1, N_HALO_FLOATS, 1, 1};

    MPI_Datatype type[9] = {MPI_INT64_T, MPI_INT, MPI_INT64_T, MPI_INT, MPI_INT,\
                            MPI_FLOAT, MPI_FLOAT, MPI_INT64_T, MPI_INT64_T};

    MPI_Type_create_struct(9, blocklen, disp, type, &this->halo_properties_MPI_Type);
    MPI_Type_commit(&this->halo_properties_MPI_Type);
  }

  {
    pid_tag_t pt;
    MPI_Aint base;
    MPI_Aint disp[3];

    MPI_Get_address(&pt, &base);
    MPI_Get_address(&pt.pid,  &disp[0]);
    MPI_Get_address(&pt.tag,  &disp[1]);
    MPI_Get_address(&pt.rank, &disp[2]);

    disp[0]-=base; disp[1]-=base; disp[2]-=base;

    int blocklen[3] = {1, 1, 1};
    MPI_Datatype type[3] = {MPI_INT64_T, MPI_INT64_T, MPI_INT};
    MPI_Type_create_struct(3, blocklen, disp, type, &this->pid_tag_MPI_Type);
    MPI_Type_commit(&this->pid_tag_MPI_Type);
  }

  // Set subextents on halo locations for this rank
  float boxStep[DIMENSION];
  for (int dim = 0; dim < DIMENSION; dim++) {
    boxStep[dim] = this->boxSize / this->layoutSize[dim];

    // define the alive region
    this->minAlive[dim] = this->layoutPos[dim] * boxStep[dim];
    this->maxAlive[dim] = this->minAlive[dim] + boxStep[dim];
    if (this->maxAlive[dim] > this->boxSize)
      this->maxAlive[dim] = this->boxSize;
  }

  if (this->myProc == MASTER) {
    cout << endl << "----------------------------------------------------------------" << endl;
    cout << "pad_factor:  " << this->pad_factor << endl;
    cout << "halo_loss_threshold:  " << this->halo_loss_threshold << endl;
    cout << "boxSize:  " << this->boxSize << endl;
    if (this->redistribute)
      cout << "REDISTRIBUTING INPUT HALO DATA" << endl;
    else
      cout << "NOT REDISTRIBUTING INPUT HALO DATA" << endl;
  }

}

void HaloDistribute::finalize() {  
  MPI_Type_free(&this->halo_properties_MPI_Type);
  MPI_Type_free(&this->pid_tag_MPI_Type);
}

void HaloDistribute::read_treenodes(string file_name) {
  
  GenericIO GIO(PartitionPlus::getComm(), file_name, GenericIO::FileIOMPI);

  // this is for restarting from a previous run -- so redistributing here is not allowed
  GIO.openAndReadHeader(GenericIO::MismatchDisallowed);
  size_t num_elem = GIO.readNumElems();

  this->halos->Resize(num_elem + GIO.requestedExtraSpace());

  GIO.addVariable("fof_halo_tag",   *(this->halos->fof_halo_tag), true);
  GIO.addVariable("fof_halo_count", *(this->halos->fof_halo_count), true);
  
  if (this->halos->has_sod)
    GIO.addVariable("sod_halo_count", *(this->halos->sod_halo_count), true);

  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names[i], *(this->halos->float_data[i]), true);

  // masses for halos and fragment halos
  GIO.addVariable("tree_node_mass", *(this->halos->tree_node_mass), true);

  // read the tree and descendant node index to populate the descendant and rank fields
  GIO.addVariable("tree_node_index", *(this->halos->tree_node_index), true);
  GIO.addVariable("desc_node_index", *(this->halos->desc_node_index), true);

  GIO.readData();

  this->halos->Resize(num_elem);

  // update the descendant index and descendant rank information
  // also validate the tree_node_index is consistent
  // set the timestep for the file -- embeded in the tree_node_index
  if (num_elem > 0) {
    for (int i=0; i<num_elem; ++i) {

      int tmp1, tmp2;
      assert(dis_tree_node_index(this->halos->tree_node_index->at(i),\
          tmp1, tmp2, this->halos->step_number));
      assert(tmp1 == i);
      assert(tmp2 == this->myProc);

      // default descendant ts to current. If there are real descendants, update.
      int ts = this->halos->step_number; 
      if (dis_tree_node_index(this->halos->desc_node_index->at(i),\
          this->halos->descendant->at(i), this->halos->rank->at(i), ts))
        this->halos->descendant_step_number = ts;
    }
  }

}

void HaloDistribute::read_halos(string file_name) {

  GenericIO GIO(PartitionPlus::getComm(), file_name, GenericIO::FileIOMPI);

  GIO.openAndReadHeader(this->MB);
  size_t num_elem = GIO.readNumElems();

  this->halos->Resize(num_elem + GIO.requestedExtraSpace());

  GIO.addVariable("fof_halo_tag",   *(this->halos->fof_halo_tag),   true);
  GIO.addVariable("fof_halo_count", *(this->halos->fof_halo_count), true);
  
  if (this->halos->has_sod)
    GIO.addVariable("sod_halo_count", *(this->halos->sod_halo_count), true);

  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names[i], *(this->halos->float_data[i]), true);

  GIO.readData();

  this->halos->Resize(num_elem);

  if (num_elem > 0) {
    // set the destination rank for the halo to unassigned
    fill(this->halos->rank->begin(), this->halos->rank->end(), this->myProc);

    // by default no descendant
    fill(this->halos->descendant->begin(), this->halos->descendant->end(), -1);

    // TODO: should I do this correctly with the cosmo params? this seems safer since this is the data
    this->halos->particle_mass = this->halos->float_data[fof_halo_mass]->at(0)/this->halos->fof_halo_count->at(0);

    // initialize the tree_node_mass to the fof_halo_mass
    copy(this->halos->float_data[fof_halo_mass]->begin(), this->halos->float_data[fof_halo_mass]->end(),\
        this->halos->tree_node_mass->begin());
  }

}

void HaloDistribute::write_halos(string file_name) {

  GenericIO GIO(PartitionPlus::getComm(), file_name);

  GIO.setPhysOrigin(0.0f);
  GIO.setPhysScale(this->boxSize);
  GIO.setNumElems(this->halos->num_halos);
  
  GIO.addVariable("fof_halo_tag", *(this->halos->fof_halo_tag));
  GIO.addVariable("fof_halo_count", *(this->halos->fof_halo_count));
  
  if (this->halos->has_sod)
    GIO.addVariable("sod_halo_count", *(this->halos->sod_halo_count));

  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names[i], *(this->halos->float_data[i]));

  GIO.write();
}

void HaloDistribute::write_pidtags(string file_name) {

  GenericIO GIO(PartitionPlus::getComm(), file_name);
  GIO.setNumElems(this->halos->pid->size());

  assert(this->halos->pid->size()==this->halos->tag->size());

  GIO.addVariable("id", this->halos->pid->begin());
  GIO.addVariable("fof_halo_tag", this->halos->tag->begin());

  GIO.write();
}

void HaloDistribute::write_treenodes(string file_name) {

  GenericIO GIO(PartitionPlus::getComm(), file_name);

  GIO.setPhysOrigin(0.0f);
  GIO.setPhysScale(this->boxSize);
  GIO.setNumElems(this->halos->num_halos);
  
  // Halo finder properties
  GIO.addVariable("fof_halo_tag", *(this->halos->fof_halo_tag));
  GIO.addVariable("fof_halo_count", *(this->halos->fof_halo_count));
  
  if (this->halos->has_sod)
    GIO.addVariable("sod_halo_count", *(this->halos->sod_halo_count));

  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names[i], *(this->halos->float_data[i]));

  // merger tree specific
  vector<int64_t> tree_node_index(this->halos->num_halos, -1);
  vector<int64_t> desc_node_index(this->halos->num_halos, -1);

  for (int i=0; i<this->halos->num_halos; ++i) {
    tree_node_index[i] = gen_tree_node_index(i, this->myProc, this->halos->step_number); 
    desc_node_index[i] = gen_tree_node_index(this->halos->descendant->at(i),\
        this->halos->rank->at(i), this->halos->descendant_step_number);
  }

  // masses for halos and fragment halos
  GIO.addVariable("tree_node_mass", *(this->halos->tree_node_mass));

  GIO.addVariable("tree_node_index", tree_node_index);
  GIO.addVariable("desc_node_index", desc_node_index);

  GIO.write();
}

void HaloDistribute::calc_halo_dest() {

  // make sure the position is inside the periodic box
  periodic_boundaries(this->halos->float_data[fof_halo_center_x],
                 this->halos->float_data[fof_halo_center_y],
                 this->halos->float_data[fof_halo_center_z],
                 this->boxSize);

  for (int i=0; i<this->halos->num_halos; ++i) {
    int home = this->myProc;
    if (this->redistribute) {
      // Figure out where to send this halo
      float xCoord = this->halos->float_data[fof_halo_center_x]->at(i);
      float yCoord = this->halos->float_data[fof_halo_center_y]->at(i);
      float zCoord = this->halos->float_data[fof_halo_center_z]->at(i);
      double sizeX = this->boxSize / this->layoutSize[0];
      double sizeY = this->boxSize / this->layoutSize[1];
      double sizeZ = this->boxSize / this->layoutSize[2];
      int coords[3] = { (int)(xCoord/sizeX),
                        (int)(yCoord/sizeY),
                        (int)(zCoord/sizeZ)
                      };

      MPI_Cart_rank(PartitionPlus::getComm(), coords, &home);

      if (home >= this->numProc) {
        ostringstream err_msg;
        err_msg << "Invalid home:" << home << " on proc " << this->myProc;
        throw std::runtime_error(err_msg.str());
      }

    } // done finding destination rank for redistribution

    this->halos->rank->at(i) = home;
  } // end of loop over every halo
}

void HaloDistribute::send_recv_halos() {
  std::vector<halo_properties_t>  haloBuffer, 
                                  haloRecvBuffer;

  std::vector<int>  haloCounts(numProc, 0), 
                    haloRecvCounts(numProc),
                    haloDisp(numProc),   
                    haloRecvDisp(numProc);

  for (int i=0; i<this->halos->num_halos; ++i) {
    ++haloCounts[this->halos->rank->at(i)];
    haloBuffer.push_back(this->halos->GetProperties(i));
  }

  this->halos->Resize(0);

  if (this->redistribute) {

    MPI_Alltoall(&haloCounts[0], 1, MPI_INT,
                 &haloRecvCounts[0], 1,MPI_INT,
                 PartitionPlus::getComm());

    size_t totalToRecv = 0;
    for (int i = 0; i < numProc; ++i)
      totalToRecv += haloRecvCounts[i];

    // Calculate receive displacements
    haloDisp[0] = haloRecvDisp[0] = 0;
    for (int i = 1; i < numProc; ++i) {
      haloDisp[i] = haloDisp[i-1] + haloCounts[i-1];
      haloRecvDisp[i] = haloRecvDisp[i-1] + haloRecvCounts[i-1];
    }

    sort(haloBuffer.begin(), haloBuffer.end(), sort_halos_by_dest);

    // Send all halos to their new homes
    haloRecvBuffer.resize(totalToRecv);
    MPI_Alltoallv(&haloBuffer[0], &haloCounts[0], &haloDisp[0], this->halo_properties_MPI_Type,
                  &haloRecvBuffer[0], &haloRecvCounts[0], &haloRecvDisp[0], this->halo_properties_MPI_Type,
                  PartitionPlus::getComm());
   
  }

  for (size_t i=0; i<haloRecvBuffer.size(); ++i)
    this->halos->PushBack(haloRecvBuffer[i]); // push back each struct variable onto corresponding vector
}

void HaloDistribute::read_pidtag(string file_name) {

  GenericIO GIO(PartitionPlus::getComm(), file_name, GenericIO::FileIOMPI);
  GIO.openAndReadHeader(this->MB); // MismatchBehavior

  size_t num_elem = GIO.readNumElems();
  this->halos->pid->resize(num_elem + GIO.requestedExtraSpace());
  this->halos->tag->resize(num_elem + GIO.requestedExtraSpace());

  // GIO treats as plain old C++ arrays
  GIO.addVariable("id"          , this->halos->pid->begin(), true);
  GIO.addVariable("fof_halo_tag", this->halos->tag->begin(), true);
  GIO.readData();

  this->halos->pid->resize(num_elem);
  this->halos->tag->resize(num_elem);

  // (b) stats about what was read by GIO
  int64_t cnt_tmp = num_elem;
  int64_t max_num_elem, min_num_elem, sum_num_elem;
  MPI_Reduce(&cnt_tmp, &max_num_elem, 1, MPI_INT64_T, MPI_MAX, 0, PartitionPlus::getComm());
  MPI_Reduce(&cnt_tmp, &min_num_elem, 1, MPI_INT64_T, MPI_MIN, 0, PartitionPlus::getComm());
  MPI_Reduce(&cnt_tmp, &sum_num_elem, 1, MPI_INT64_T, MPI_SUM, 0, PartitionPlus::getComm());

  if (this->myProc == MASTER) {
    cout << "Local PID-TAG read stats" << endl;
    cout << "max:" << max_num_elem << " min:" << min_num_elem << " avg:" 
      << sum_num_elem/PartitionPlus::getNumProc() << endl;
  }
}

void HaloDistribute::send_recv_pidtag() {

  // pack on the second half of the buffer
  allocated_vector<pid_tag_t> send_buff;
  send_buff.set_buffer((void*)(this->halos->buffer + this->halos->buffer_size/2));
  send_buff.set_capacity(this->max_num_pidtags);

  for ( size_t i = 0; i < this->halos->pid->size(); ++i ) {
    pid_tag_t pt = { this->halos->pid->at(i), this->halos->tag->at(i),\
        this->halos->rank->at((*this->halos->tag2idx)[this->halos->tag->at(i)]) };
    send_buff.push_back(pt);
  }

  // overwrite the GIO read buffers on the first half
  memcpy((void*)this->halos->buffer, (const void*)send_buff.begin(), send_buff.size()*sizeof(pid_tag_t));
  send_buff.set_buffer(this->halos->buffer);

  // MPI receive on the second half of the buffer
  allocated_vector<pid_tag_t> recv_buff;
  recv_buff.set_buffer((void*)(this->halos->buffer + this->halos->buffer_size/2));
  recv_buff.set_capacity(this->max_num_pidtags);

  vector<int> pidtag_send_counts(this->numProc,0);
  for (int i=0; i<this->halos->num_halos; ++i)
    pidtag_send_counts.at(this->halos->rank->at(i)) += this->halos->fof_halo_count->at(i);

  if (this->redistribute) {
    
    sort(send_buff.begin(), send_buff.end(), sort_by_dest);

    for (size_t i=0; i<send_buff.size(); ++i) {
      if (send_buff.at(i).rank>= PartitionPlus::getNumProc())
        cerr << "Destination error! " << send_buff.at(i).rank << " on proc " << this->myProc << endl;
    }

    std::vector<int> hptRecvCounts(numProc, 0),
                     hptDisp(numProc, 0),
                     hptRecvDisp(numProc, 0);

    // Record all sizes into a single buffer and send this to all ranks
    size_t totalToSend = 0;
    for (int i = 0; i < this->numProc; ++i)
      totalToSend += pidtag_send_counts[i];

    MPI_Alltoall( &pidtag_send_counts[0], 1, MPI_INT,
                  &hptRecvCounts[0], 1, MPI_INT,
                  PartitionPlus::getComm());

    // get from every other rank
    size_t totalToRecv = 0;
    for (int i = 0; i < numProc; ++i) {
      totalToRecv += hptRecvCounts[i];
    }

    // calculate displacements
    hptDisp[0] = hptRecvDisp[0] = 0;
    for (int i = 1; i < numProc; ++i) {
      hptDisp[i] = hptDisp[i-1] + pidtag_send_counts[i-1];
      hptRecvDisp[i] = hptRecvDisp[i-1] + hptRecvCounts[i-1];
    }

    // Send all (pid,tag)s to their new homes
    recv_buff.resize(totalToRecv);

    int max_send = 0;
    std::vector<int> hptCounts_t (numProc, 0),
                     hptRecvCounts_t(numProc, 0);

    for (int i = 0; i < numProc; ++i) {
      if (max_send < pidtag_send_counts[i])
        max_send = pidtag_send_counts[i];
    }

    int max_send_all;
    MPI_Allreduce(&max_send, &max_send_all, 1, MPI_INT, MPI_MAX, PartitionPlus::getComm());
    if (this->myProc == MASTER)
      cout << "Maximum rank-to-rank message size: " 
        << max_send_all*sizeof(pid_tag_t)/GB_BYTES << " GB -- limit is " 
        << MSG_SIZE/GB_BYTES << " GB" << endl;
      
    const int max_send_items = MSG_SIZE/sizeof(pid_tag_t);
    int numParts = max_send_all/max_send_items + (int)((max_send_all % max_send_items) > 0);
    if (this->myProc == MASTER)
      cout << "Sending in "  << numParts << " iterations." << endl;
    
    for (int i = 0; i < numParts; ++i) {
      for (int j = 0; j < numProc; ++j) {
        hptCounts_t[j]     = (pidtag_send_counts[j] > max_send_items) ? max_send_items : pidtag_send_counts[j];
        hptRecvCounts_t[j] = (hptRecvCounts[j] > max_send_items) ? max_send_items : hptRecvCounts[j];
      }
      MPI_Alltoallv( &send_buff[0], &hptCounts_t[0], &hptDisp[0], this->pid_tag_MPI_Type,
                     &recv_buff[0], &hptRecvCounts_t[0], &hptRecvDisp[0], this->pid_tag_MPI_Type,
                     PartitionPlus::getComm());
      for (int j = 0; j < numProc; ++j) {
        hptDisp[j]     += hptCounts_t[j];
        hptRecvDisp[j] += hptRecvCounts_t[j];
        pidtag_send_counts[j] -= hptCounts_t[j];
        hptRecvCounts[j]  -= hptRecvCounts_t[j];
      }
    }

  }
  else { // no redistribute
    // just copy the send_buff to the recv_buff
    memcpy((void*)recv_buff.begin(), (const void*)send_buff.begin(), send_buff.size()*sizeof(pid_tag_t));
  }

  // We now have all of our (pid,tag) pairs, put them in our local arrays
  this->halos->pid->set_capacity(recv_buff.size());
  this->halos->pid->resize(0);

  this->halos->tag->set_capacity(recv_buff.size());
  this->halos->tag->resize(0);

  for (size_t i = 0; i < recv_buff.size(); ++i) {
    this->halos->pid->push_back(recv_buff[i].pid);
    this->halos->tag->push_back(recv_buff[i].tag);
  }

}

void HaloDistribute::calc_max_num_pidtag() {
  
  this->max_num_pidtags = 0;

  // count the number of pidtag pairs sending to each rank
  vector<int> pidtag_send_counts(this->numProc,0);
  vector<int> pidtag_recv_counts(this->numProc,0);

  for (int i=0; i<this->halos->num_halos; ++i) {
    pidtag_send_counts.at(this->halos->rank->at(i)) += this->halos->fof_halo_count->at(i);

    // count what will be read locally
    this->max_num_pidtags += this->halos->fof_halo_count->at(i);
  }

  MPI_Allreduce(&pidtag_send_counts[0], &pidtag_recv_counts[0], this->numProc, MPI_INT, MPI_SUM, PartitionPlus::getComm()); 

  // now we know the number of PIDTAG pairs that will be on each rank
  if (this->max_num_pidtags < pidtag_recv_counts[this->myProc])
    this->max_num_pidtags = pidtag_recv_counts[this->myProc];

  // ** duplicate removal with neghbors uses the same memory buffers
  const size_t send_cap = MSG_SIZE/sizeof(pid_tag_t);

  for (int i=0; i<NUM_OF_NEIGHBORS; ++i) {

    const int recv_cnt_ngh = pidtag_recv_counts[this->neighbor[i]];

    if (send_cap > recv_cnt_ngh && recv_cnt_ngh > this->max_num_pidtags)
      this->max_num_pidtags = recv_cnt_ngh;

    else if (send_cap < recv_cnt_ngh && send_cap > this->max_num_pidtags)
      this->max_num_pidtags = send_cap;
  }

}

void HaloDistribute::allocate_pidtag_buffer() {
  if (this->halos->buffer==NULL) {

    this->halos->buffer_size = 2 * (size_t)(this->pad_factor * this->max_num_pidtags * sizeof(pid_tag_t));
    this->halos->buffer = (char*)malloc(this->halos->buffer_size);

    if (this->halos->buffer==NULL) {
      ostringstream err_msg;
      err_msg << "rank:" << this->myProc << " could not allocate " << 
        this->halos->buffer_size/GB_BYTES << "GB of memory." << endl;
      throw std::runtime_error(err_msg.str());
    }
  }
  else { // buffer is pre-allocated
    if (this->max_num_pidtags > this->halos->buffer_size/sizeof(pid_tag_t)/2.f) {
      ostringstream err_msg;
      err_msg << "Insufficient buffer size!" << endl;
      throw std::runtime_error(err_msg.str());
    }
  }

  // point the pid and tag allocated vectors to the new buffer
  this->halos->pid->set_buffer((void*)this->halos->buffer);
  this->halos->pid->set_capacity(this->max_num_pidtags * this->pad_factor);
  this->halos->tag->set_buffer((void*)this->halos->pid->cap_end());
  this->halos->tag->set_capacity(this->max_num_pidtags * this->pad_factor);
}

void HaloDistribute::remove_duplicate_PIDs() {

  // convert the parallel arrays: PID and TAG to an array of structures

  size_t removed_local = 0;

  allocated_vector<pid_tag_t> send_buff;
  send_buff.set_buffer(this->halos->buffer + this->halos->buffer_size/2);
  send_buff.set_capacity(this->max_num_pidtags);

  if (this->halos->pid->size() > 0) {

    // NOTE: We're abusing this structure to hold the fof halo count (in the rank field),
    // which is use to make a decision on which particle to delete.
    for (size_t i = 0; i < this->halos->pid->size(); ++i) {
      pid_tag_t pt = { this->halos->pid->at(i), this->halos->tag->at(i),\
          this->halos->fof_halo_count->at((*this->halos->tag2idx)[this->halos->tag->at(i)]) };
      send_buff.push_back(pt);
    }

    this->halos->pid->resize(0);
    this->halos->tag->resize(0);

    // overwrite the GIO read buffers on the first half
    memcpy((void*)this->halos->buffer, (const void*)send_buff.begin(), send_buff.size()*sizeof(pid_tag_t));
    send_buff.set_buffer(this->halos->buffer);

    // the rest of this function assumes PIDs are sorted
    sort(send_buff.begin(), send_buff.end(), sort_by_pid);

    // 
    // *** remove local duplicates (on the same rank) ***
    //

    for (size_t i = 0; i < send_buff.size()-1; ++i) {
      if (send_buff[i].pid == send_buff[i+1].pid) {
        // NOTE: destination is actually the fof halo count here!
        // ^^^^  **********************************************
        // if the halo counts are not equal -- rank field is holding count
        if (send_buff[i].rank != send_buff[i+1].rank) {
          // if right side is bigger, keep the left
          //if (send_buff[i].rank < send_buff[i+1].rank) {
          if (send_buff[i].rank > send_buff[i+1].rank) {
            swap(send_buff[i+1], send_buff[i]);  // overwrite the one on the right
          }
        }
        else { // the halo counts are equal
          // we need to choose -- and it needs to be consistent
          // we'll keep the one with a lower halo tag
          // NOTE: if tags are equal we already ensure that halos on the
          // same rank have unique tag so it dosen't matter which one we keep
          if (send_buff[i].tag < send_buff[i+1].tag) {
            swap(send_buff[i+1], send_buff[i]);
          }
        }
        send_buff[i].pid = -1;
        ++removed_local;
      }
    }

  } // end if > 0 pids

  // get the total count from all ranks
  size_t removed_total;
  MPI_Reduce(&removed_local, &removed_total, 1, MPI_INT64_T, MPI_SUM, MASTER, PartitionPlus::getComm());
  if (this->myProc == MASTER)
    cout << "Total PIDs locally removed:" << removed_total << endl;

  //
  // *** remove duplicates on neighboring ranks ***
  //

  // MPI receive on the second half of the buffer
  allocated_vector<pid_tag_t> recv_buff;
  recv_buff.set_buffer((void*)(this->halos->buffer + this->halos->buffer_size/2));
  recv_buff.set_capacity(this->max_num_pidtags);

  size_t block_size;
  size_t cap = min((size_t)(MSG_SIZE/sizeof(pid_tag_t)), recv_buff.capacity());

  MPI_Allreduce(&cap, &block_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, PartitionPlus::getComm());

  size_t removed_remote = 0;
  removed_total = 0;
  // loop over the neighbor ranks
  for (int i=0; i<NUM_OF_NEIGHBORS; i++)
    removed_remote += this->neighbor_exchange(i, send_buff, recv_buff, block_size);

  MPI_Reduce(&removed_remote, &removed_total, 1, MPI_INT64_T, MPI_SUM, MASTER, PartitionPlus::getComm());
  if (this->myProc == MASTER)
    cout << "Total PIDs remotely removed:" << removed_total << endl;

  //
  // *** update the halos ***
  //

  // copy the original fof_halo_count -- used to compare 
  vector<int> original_fof_halo_count = (*this->halos->fof_halo_count);

  // adjust the halo counts for the PIDTAGs we remove
  for (size_t i = 0; i < send_buff.size(); ++i) {
    if (send_buff[i].pid==-1)
      --this->halos->fof_halo_count->at((*this->halos->tag2idx)[send_buff[i].tag]);
  }

  // filter out any halos that are now below the threshold
  // particles in removed halos are flagged by (-1) as the duplicates
  size_t local_removed = this->remove_halos_below_threshold(send_buff, original_fof_halo_count);
  size_t total_removed;
  MPI_Reduce(&local_removed, &total_removed, 1, MPI_INT64_T, MPI_SUM, MASTER, PartitionPlus::getComm());
  if (this->myProc==MASTER)
    cout << "Removed " << total_removed << " halos that lost more than " << 100-this->halo_loss_threshold*100 << "\% of their particles." << endl;

  //
  // *** actually remove flagged particles ***
  //

  if (send_buff.size() > 0) {

    // move all the marked (-1) duplicates to the begining of the vector
    sort(send_buff.begin(), send_buff.end(), sort_by_pid);

    size_t removed_all = 0;
    for (size_t i = 0; i < send_buff.size(); ++i) {
      if (send_buff[i].pid==-1)
        ++removed_all;
    }

    // chop off the first few we are removing
    send_buff.set_buffer((void*)(send_buff.begin() + removed_all));
    send_buff.resize(send_buff.size() - removed_all);

    //
    // *** unpack the PIDTAGs from SOA to parallel arraysf ***
    //

    // copy the send buffer to the address at the end because pid and hid are at the begining
    // memcpy(dest, src, size)
    memcpy((void*)(this->halos->buffer + this->halos->buffer_size/2), (const void*)send_buff.begin(), send_buff.size()*sizeof(pid_tag_t));
    send_buff.set_buffer((void*)(this->halos->buffer + this->halos->buffer_size/2));

    for (int i=0; i<send_buff.size(); ++i) {
      this->halos->pid->push_back(send_buff[i].pid);
      this->halos->tag->push_back(send_buff[i].tag);
    }

  } // end if send_buff.size() > 0

  // 
  // *** update the FOF halo masses ***
  //

  this->halos->RecalculateFoFmass();

}

size_t HaloDistribute::remove_halos_below_threshold(allocated_vector<pid_tag_t> &pidtag_buff, vector<int> original_fof_halo_count) {

  size_t local_removed = this->halos->num_halos;

  // pack into a AOS buffer 
  std::vector<halo_properties_t> haloBuffer;

  for (int i=0; i<this->halos->num_halos; ++i) {
    //if (this->halos->fof_halo_count->at(i) >= this->min_halo_count) // don't pack below threshold halos
    if (this->halos->fof_halo_count->at(i) >= original_fof_halo_count.at(i)*this->halo_loss_threshold) // don't pack below threshold halos
      haloBuffer.push_back(this->halos->GetProperties(i));
    //else
    //  cout << "Deleting halo with " << original_fof_halo_count.at(i) << "-->" << this->halos->fof_halo_count->at(i) << " particles." << endl;
  }

  // unpack to SOA
  this->halos->Resize(0);
  for (size_t i=0; i<haloBuffer.size(); ++i)
    this->halos->PushBack(haloBuffer[i]); // push back each struct variable onto corresponding vector

  this->halos->VerifyHalos();

  // some particles may have become orphaned -- delete them if not in the new halo map
  for (int i=0; i<pidtag_buff.size(); ++i)
    if (this->halos->tag2idx->count(pidtag_buff[i].tag)==0)
      pidtag_buff[i].pid=-1;

  local_removed -= this->halos->num_halos;

  return local_removed;
}

// on big memory machines the amount of data that needs to be sent to other ranks
// can exceed 32-bit MPI counts or make the buffers prohibitively large, so 
// a block size (in bytes) can be specified to drive an iterative approach. 
int HaloDistribute::neighbor_exchange(int neighbor,\
    allocated_vector<pid_tag_t> &send_buff, allocated_vector<pid_tag_t> &recv_buff,\
    size_t block_size) {
    
  int num_deleted = 0;

  int send_idx = neighbor;
  int recv_idx = neighbor % 2 ? send_idx - 1 : send_idx + 1;

  MPI_Request send_request;
  MPI_Status  status;

  // loop until there is nothing to send or receive
  bool sending = true;
  bool receiving = true;
  size_t send_offset = 0;
  size_t send_size = 0;
  int num_recv;

  while (sending || receiving) {
    size_t remaining = 0;
    if (send_offset <= send_buff.size())
      remaining = send_buff.size() - send_offset;

    send_size = sending ? (remaining < block_size ? remaining : block_size) : 0;
    if (send_size==0)
      send_offset = 0;

    if (sending) {
      pid_tag_t *ptr = NULL; 
      if (send_size > 0)
        ptr = &send_buff.at(send_offset);
      MPI_Isend(ptr, send_size, pid_tag_MPI_Type, this->neighbor[send_idx], 0, PartitionPlus::getComm(), &send_request );
    }

    if (receiving) {
      MPI_Probe(this->neighbor[recv_idx], 0, PartitionPlus::getComm(), &status);
      MPI_Get_count(&status, pid_tag_MPI_Type, &num_recv);
      recv_buff.resize(num_recv);
      if (num_recv==0)
        MPI_Recv(NULL, num_recv, pid_tag_MPI_Type, this->neighbor[recv_idx], 0, PartitionPlus::getComm(), MPI_STATUS_IGNORE);
      else {

        MPI_Recv(&recv_buff.at(0), num_recv, pid_tag_MPI_Type, this->neighbor[recv_idx], 0, PartitionPlus::getComm(), MPI_STATUS_IGNORE);

        // we can assume the PIDs are unique on each rank here
        int i, j;
        i = j = 0;
        while (i<send_buff.size() && j<recv_buff.size()) {
          if (send_buff[i].pid==-1)
            ++i;
          else if (recv_buff[j].pid==-1)
            ++j;
          else if (send_buff[i].pid==recv_buff[j].pid) {
            // here I am abusing the rank field again to hold the halo count
            if (send_buff[i].rank != recv_buff[j].rank) {   // not equal halo count
              //if (send_buff[i].rank > recv_buff[j].rank) {  // keep the particle on the SMALLER halo
              if (send_buff[i].rank < recv_buff[j].rank) {  // keep the particle on the LARGER halo
                send_buff[i].pid=-1;
                ++num_deleted;
              }
            }
            else { // same PID equal halo count
              if (send_buff[i].tag != recv_buff[j].tag) {   // not equal halo tag
                if (send_buff[i].tag > recv_buff[j].tag) {  // keep the particle with a smaller tag
                  send_buff[i].pid=-1;
                  ++num_deleted;
                }
              }
              else { // equal halo count and equal tag
                if (this->myProc > this->neighbor[recv_idx]) {  // keep the particle from a lower rank
                  send_buff[i].pid=-1;
                  ++num_deleted;
                }
              }
            }
            ++i;
            ++j;
          }
          else if (send_buff[i].pid<recv_buff[j].pid)
            ++i;
          else if (send_buff[i].pid>recv_buff[j].pid)
            ++j;
        }
      }
    }

    if (sending)
      MPI_Wait( &send_request, MPI_STATUS_IGNORE );

    send_offset += send_size;
    sending = send_size > 0;
    receiving = num_recv > 0; 
  } // end of while loop for sending blocks

  return num_deleted;
}

void HaloDistribute::read_treenodes_pidtags_restart(string baseFileTreenodes,\
    string baseFileHaloParticleTags, int step_number) {

  PartitionPlus::programStatus("Reading treenodes from step.");
  this->halos->step_number = step_number;
  this->halos->descendant_step_number = step_number;

  PartitionPlus::programStatus("Reading treenodes from GIO.");
  this->read_treenodes(baseFileTreenodes);

  PartitionPlus::programStatus("Verifying halos");
  this->halos->VerifyHalos();

  PartitionPlus::programStatus("Calculating PIDTAG memory buffer size");
  this->calc_max_num_pidtag(); 

  PartitionPlus::programStatus("Allocating PIDTAG memory buffer");
  this->allocate_pidtag_buffer();

  PartitionPlus::programStatus("Reading PIDTAGs from GIO");
  this->read_pidtag(baseFileHaloParticleTags);

  PartitionPlus::programStatus("Verifying PIDTAGs");
  this->halos->VerifyPIDTAGs();

}

// does the GIO read and redistribution of halos and their particle IDs
void HaloDistribute::read_halos_pidtags(const string baseFileFofProperties,\
    const string baseFileHaloParticleTags, const int step_number) {

  // TIMERS
  SimpleTimings::TimerRef t_all = SimpleTimings::getTimer("ALL_timer");
  SimpleTimings::TimerRef t_gio = SimpleTimings::getTimer("GIO_timer");
  SimpleTimings::TimerRef t_net = SimpleTimings::getTimer("NET_timer");

  SimpleTimings::startTimer(t_all);

  char step[8];
  int n = sprintf(step, "%d", step_number);
  PartitionPlus::programStatus("Reading halos from step " + string(step));

  // set the step number 
  this->halos->step_number = step_number;
  this->halos->descendant_step_number = step_number;

  PartitionPlus::programStatus("Reading halo properties...");
  SimpleTimings::startTimer(t_gio);
    // Read halo properties from GIO
    this->read_halos(baseFileFofProperties);
  SimpleTimings::stopTimer(t_gio);

  // using the halo position, compute the ranks halos should live on
  // store the destination rank this->halos->at(i).rank
  this->calc_halo_dest();

  // calc the max to read OR receive from another process
  // needs calc_halo_dest to be called to know the redistribution
  this->calc_max_num_pidtag(); 

  PartitionPlus::programStatus("Allocating PIDTAG buffer...");
  this->allocate_pidtag_buffer();

  PartitionPlus::programStatus("Reading PIDTAGs...");
  SimpleTimings::startTimer(t_gio);
    this->read_pidtag(baseFileHaloParticleTags);
  SimpleTimings::stopTimer(t_gio);

  // READING FILES COMPLETE

  size_t tot_n_halo;

  MPI_Reduce(&this->halos->num_halos, &tot_n_halo,\
      1, MPI_INT64_T, MPI_SUM, cosmotk::MASTER,\
      cosmotk::PartitionPlus::getComm());

  char n_halo[32];
  char n_dup[32];

  string message;

  sprintf(n_halo, "%ld", tot_n_halo);
  // in case of duplicate halos, only the first halo is kept
  sprintf(n_dup, "%ld", this->halos->VerifyHalos()); // return value is dup count from all ranks
  PartitionPlus::programStatus(string(n_dup) + " duplicate fof_halo_tag values out of " + string(n_halo) + " total halos");

  sprintf(n_dup, "%ld", this->halos->VerifyPIDTAGs()); 
  message = "All particles belong to a halo.\nAll halos have the correct number of particles.\n" + string(n_dup) + " particles were deleted.";
  PartitionPlus::programStatus(message);

  // this just sends the PIDs to the rank where the tag matches locally
  // ... uses the result from calc_halo_dest
  PartitionPlus::programStatus("Distributing PIDTAGs...");
  this->send_recv_pidtag();

  PartitionPlus::programStatus("Distributing halos...");
  this->send_recv_halos();

  MPI_Reduce(&this->halos->num_halos, &tot_n_halo,\
      1, MPI_INT64_T, MPI_SUM, cosmotk::MASTER,\
      cosmotk::PartitionPlus::getComm());

  sprintf(n_halo, "%ld", tot_n_halo);
  // in case of duplicate halos, only the first halo is kept
  sprintf(n_dup, "%ld", this->halos->VerifyHalos()); // return value is dup count from all ranks
  PartitionPlus::programStatus(string(n_dup) + " duplicate fof_halo_tag values out of " + string(n_halo) + " total halos");

  sprintf(n_dup, "%ld", this->halos->VerifyPIDTAGs()); 
  message = "All particles belong to a halo.\nAll halos have the correct number of particles.\n" + string(n_dup) + " particles were deleted.";
  PartitionPlus::programStatus(message);

  /*PartitionPlus::programStatus("Reverifying halos...");
  this->halos->VerifyHalos();

  PartitionPlus::programStatus("Reverifying PIDTAGs...");
  this->halos->VerifyPIDTAGs();*/

  PartitionPlus::programStatus("Removing duplicate PIDs and any halos that fall below threshold.");
  this->remove_duplicate_PIDs();

  // tree_node_mass is initially fof_halo_mass, but can become different due to fragments
  // copy fof_halo_mass to tree_node_mass
  std::copy(this->halos->float_data[fof_halo_mass]->begin(), this->halos->float_data[fof_halo_mass]->end(),\
      this->halos->tree_node_mass->begin());

  PartitionPlus::programStatus("Reverifying PIDTAGs after duplicate removal.");
  this->halos->VerifyPIDTAGs();

  SimpleTimings::stopTimer(t_all);

  SimpleTimings::timerStats(t_net);
  SimpleTimings::timerStats(t_gio);
  SimpleTimings::timerStats(t_all);

  nodeMemory(string("final").c_str(),true,GB);
}

} // END namespace cosmotk

  // test breaking the verifications
  //this->pid->resize(this->pid->size()-10);
  //this->hid.resize(this->hid.size()-10);
  //this->hid.at(99)=-200;

