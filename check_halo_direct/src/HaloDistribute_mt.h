#ifndef HALO_DISTRIBUTE_MT_H
#define HALO_DISTRIBUTE_MT_H

#include <string>
#include <stdint.h>
#include "Halos_test.h"
#include "AllocatedVector.h"

using std::string;
using std::ifstream;
using std::vector;

struct pid_tag_t {
  int64_t pid;
  int64_t tag;
  int rank;
};

namespace cosmotk {

  class HaloDistribute {

  public:
    HaloDistribute(Halos*);
    ~HaloDistribute();

    void set_parameters(                      // Set parameters halo distribution
        bool redistribute,                    // sets GIO mismatch behaviour and some program control flow
        float pad_factor,                     // memory extra space
        POSVEL_T rL,                          // Box size of the physical problem
        float halo_loss_threshold);           // halos below this threshold are filtered out

    void initialize();
    void read_halos_pidtags(const string, const string, const int);
    void read_treenodes_pidtags_restart(string, string, int);
    void write_halos(string file_name);
    void write_pidtags(string file_name);
    void read_treenodes(string);
    void write_treenodes(string file_name);
    void finalize();
    void read_halos(string);

    MPI_Datatype halo_properties_MPI_Type;    //

  private:
    gio::GenericIO::MismatchBehavior MB;      // sets if GIO will do n-to-n read or n-to-m
    bool    redistribute;                     //  
    float   boxSize;                          // Physical box size (rL)
    float   pad_factor;                       // 
    int     myProc;                           // My processor number
    int     numProc;                          // Total number of processors
    int     layoutSize[DIMENSION];            // Decomposition of processors
    int     layoutPos[DIMENSION];             // Position of this processor in decomposition
    float   halo_loss_threshold;              // filter out halos that have lost more than this threshold * particle count
    size_t  max_num_pidtags;
    float   minAlive[DIMENSION];              // Minimum alive halo location on processor
    float   maxAlive[DIMENSION];              // Maximum alive halo location on processor
    int     neighbor[NUM_OF_NEIGHBORS];       // Neighbor processor ids
    Halos* halos;                             // the halo data class
    MPI_Datatype pid_tag_MPI_Type;            //

    //void read_halos(string);
    void read_pidtag(string);
    void calc_halo_dest();
    void calc_max_num_pidtag();
    void allocate_pidtag_buffer(); 
    void send_recv_halos();
    void send_recv_pidtag();
    void remove_duplicate_PIDs();
    int neighbor_exchange(int, allocated_vector<pid_tag_t>&, allocated_vector<pid_tag_t>&, size_t);
    size_t remove_halos_below_threshold(allocated_vector<pid_tag_t>&, vector<int>);
    //size_t remove_halos_below_threshold(allocated_vector<pid_tag_t>&);
  };
}
#endif
