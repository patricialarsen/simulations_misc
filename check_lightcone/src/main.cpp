#include <cstdlib>
#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdexcept>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <algorithm>    // std::sort
#include <sstream>
#include <omp.h>
#include <sys/stat.h>

// Generic IO
#include "GenericIO.h"
#include "Partition.h"


#include "MurmurHashNeutral2.cpp" 


/// Assumes ids and replication values are consistent but ordering may have changed
/// Redistributes amongst ranks and sorts the values for a one-to-one check of the changes in 
/// each of the main lightcone outputs 


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;


typedef struct lc_halo {
  float posvel_a_m[7];
  float phi;
  int rr;
  int64_t id;
  unsigned int destination_rank;
} lc_halo;

bool comp_by_lc_dest(const lc_halo &a, const lc_halo &b) {
  return a.destination_rank < b.destination_rank;
}

bool comp_by_id(const lc_halo &a, const lc_halo &b) {
  return a.id < b.id;
}

bool comp_by_rep(const lc_halo &a, const lc_halo &b) {
  return a.rr < b.rr;
}

void  compute_mean_std_dist(vector<float> val1 , vector<float> val2 , float diff_tot ,string var_name){
  // compute the mean and std of the differences and make a histogram to look more closely
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = val1.size();
  double diff=0; 
  double diff_frac=0;
  int n_tot;
  double frac_max;
  double mean;
  double stddev=0;

  for (int i=0; i<count; i++){
      diff += (double)(val1[i]-val2[i]);
      if (val1[i]!=0){
        double frac = (double)(fabs(val1[i]-val2[i])/fabs(val1[i]));
         diff_frac = (diff_frac<frac)?frac:diff_frac;
      }
   }

      MPI_Allreduce(&diff, &mean, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
      MPI_Reduce(&diff_frac, &frac_max, 1, MPI_DOUBLE, MPI_MAX, 0, Partition::getComm());
      MPI_Allreduce(&count, &n_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());  
   mean = mean/n_tot;   

   for (int i=0; i< count; i++){
      stddev += (double) ((val1[i]-val2[i]-(float)mean)*(val1[i]-val2[i]-(float)mean)/(n_tot-1));
   }
   double stddev_tot;
   MPI_Reduce(&stddev, &stddev_tot, 1, MPI_DOUBLE, MPI_SUM, 0, Partition::getComm());
   stddev_tot = sqrt(stddev_tot);

   if (rank==0){
     cout << var_name << endl;
     cout << "______________________________________" <<endl;
     cout << " mean difference = "<< mean << endl;
     cout << " maximum fractional difference = "<< frac_max<< endl;
     cout << " standard deviation of difference = " << stddev_tot << endl;
     cout << endl;
   }



  return;
}




struct IO_Buffers { // this is the order of lc_halo struct above

  // LC halo data
  vector<float> x;
  vector<float> y;
  vector<float> z;
  vector<float> vx;
  vector<float> vy;
  vector<float> vz;
  vector<float> a;
  vector<float> phi;
  vector<int> replication;
  vector<int64_t> id;
  double box_size[3];
  double origin[3];

  // 
};


IO_Buffers IOB;
IO_Buffers IOB2;

void erase_element_IOB2(int i){
  IOB2.x.erase(IOB2.x.begin()+i);
  IOB2.y.erase(IOB2.y.begin()+i);
  IOB2.z.erase(IOB2.z.begin()+i);
  IOB2.vx.erase(IOB2.vx.begin()+i);
  IOB2.vy.erase(IOB2.vy.begin()+i);
  IOB2.vz.erase(IOB2.vz.begin()+i);
  IOB2.a.erase(IOB2.a.begin()+i);
  IOB2.replication.erase(IOB2.replication.begin()+i);
  IOB2.id.erase(IOB2.id.begin()+i);
  IOB2.phi.erase(IOB2.phi.begin()+i);
}


void erase_element_IOB(int i){
  IOB.x.erase(IOB.x.begin()+i);
  IOB.y.erase(IOB.y.begin()+i);
  IOB.z.erase(IOB.z.begin()+i);
  IOB.vx.erase(IOB.vx.begin()+i);
  IOB.vy.erase(IOB.vy.begin()+i);
  IOB.vz.erase(IOB.vz.begin()+i);
  IOB.a.erase(IOB.a.begin()+i);
  IOB.replication.erase(IOB.replication.begin()+i);
  IOB.id.erase(IOB.id.begin()+i);
  IOB.phi.erase(IOB.phi.begin()+i);
}

void clear_IO_buffers() {
  IOB.x.clear();
  IOB.y.clear();
  IOB.z.clear();
  IOB.vx.clear();
  IOB.vy.clear();
  IOB.vz.clear();
  IOB.a.clear();
  IOB.replication.clear();
  IOB.id.clear();
  IOB.phi.clear();
  IOB2.x.clear();
  IOB2.y.clear();
  IOB2.z.clear();
  IOB2.vx.clear();
  IOB2.vy.clear();
  IOB2.vz.clear();
  IOB2.a.clear();
  IOB2.replication.clear();
  IOB2.id.clear();
  IOB2.phi.clear();

}

inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}


void read_lc_file(string file_name, string file_name2) {
  // reads in relevant lightcone information from a particle or halo file 
  {
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  GIO.readPhysScale(IOB.box_size);
  GIO.readPhysOrigin(IOB.origin);
  IOB.x.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.y.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.z.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.vx.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.vy.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.vz.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.a.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB.replication.resize(num_elems + GIO.requestedExtraSpace()/sizeof(int));
  IOB.id.resize(num_elems + GIO.requestedExtraSpace()/sizeof(int64_t));
  IOB.phi.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  GIO.addVariable("x", IOB.x, true);
  GIO.addVariable("y", IOB.y, true);
  GIO.addVariable("z", IOB.z, true);
  GIO.addVariable("vx", IOB.vx, true);
  GIO.addVariable("vy", IOB.vy, true);
  GIO.addVariable("vz", IOB.vz, true);
  GIO.addVariable("a", IOB.a, true);
  GIO.addVariable("replication", IOB.replication, true);
  GIO.addVariable("id", IOB.id, true);
  GIO.addVariable("phi",IOB.phi,true);
  GIO.readData();
  IOB.x.resize(num_elems);
  IOB.y.resize(num_elems);
  IOB.z.resize(num_elems);
  IOB.vx.resize(num_elems);
  IOB.vy.resize(num_elems);
  IOB.vz.resize(num_elems);
  IOB.a.resize(num_elems);
  IOB.replication.resize(num_elems);
  IOB.id.resize(num_elems);
  IOB.phi.resize(num_elems);
  }
  {
  GenericIO GIO(Partition::getComm(),file_name2,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  GIO.readPhysScale(IOB2.box_size);
  GIO.readPhysOrigin(IOB2.origin);
  IOB2.x.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB2.y.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB2.z.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB2.vx.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB2.vy.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB2.vz.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB2.a.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));
  IOB2.replication.resize(num_elems + GIO.requestedExtraSpace()/sizeof(int));
  IOB2.id.resize(num_elems + GIO.requestedExtraSpace()/sizeof(int64_t));
  IOB2.phi.resize(num_elems + GIO.requestedExtraSpace()/sizeof(float));

  GIO.addVariable("x", IOB2.x, true);
  GIO.addVariable("y", IOB2.y, true);
  GIO.addVariable("z", IOB2.z, true);
  GIO.addVariable("vx", IOB2.vx, true);
  GIO.addVariable("vy", IOB2.vy, true);
  GIO.addVariable("vz", IOB2.vz, true);
  GIO.addVariable("a", IOB2.a, true);
  GIO.addVariable("replication", IOB2.replication, true);
  GIO.addVariable("id", IOB2.id, true);
  GIO.addVariable("phi", IOB2.phi, true);

  GIO.readData();
  IOB2.x.resize(num_elems);
  IOB2.y.resize(num_elems);
  IOB2.z.resize(num_elems);
  IOB2.vx.resize(num_elems);
  IOB2.vy.resize(num_elems);
  IOB2.vz.resize(num_elems);
  IOB2.a.resize(num_elems);
  IOB2.phi.resize(num_elems);

  IOB2.replication.resize(num_elems);
  IOB2.id.resize(num_elems);
 }
}


int main( int argc, char** argv ) {
  MPI_Init( &argc, &argv );
  Partition::initialize();
  GenericIO::setNaturalDefaultPartition();

  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  string lc_file     = string(argv[1]);
  string lc_file2    = string(argv[2]);

  int err = 0;
  bool skip_err = true;

  // LC HALOS
  read_lc_file(lc_file, lc_file2);                  // read the file

  if (rank == 0)
    cout << "Done reading LC" << endl;

  vector<lc_halo> lc_halo_send;
  vector<int> lc_halo_send_cnt;
  lc_halo_send_cnt.resize(n_ranks,0);
  for (size_t i=0;i<IOB.id.size();++i) {     // pack a vector of structures
    int64_t tag = IOB.id[i];
    unsigned int recv_rank = tag_to_rank(tag, n_ranks);

    lc_halo h = { {IOB.x[i],IOB.y[i],IOB.z[i],IOB.vx[i],IOB.vy[i],IOB.vz[i],IOB.a[i]}, \
                  IOB.phi[i],IOB.replication[i],\
                  tag, recv_rank};
    lc_halo_send.push_back(h);
    ++lc_halo_send_cnt[recv_rank];        // count the send items
  }                                       // prepare to send by contiguous destinations


  vector<lc_halo> lc_halo_send2;
  vector<int> lc_halo_send_cnt2;
  lc_halo_send_cnt2.resize(n_ranks,0);
  for (size_t i=0;i<IOB2.id.size();++i) {     // pack a vector of structures
    int64_t tag = IOB2.id[i];
    unsigned int recv_rank = tag_to_rank(tag, n_ranks);

    lc_halo h = { {IOB2.x[i],IOB2.y[i],IOB2.z[i],IOB2.vx[i],IOB2.vy[i],IOB2.vz[i],IOB2.a[i]}, \
                  IOB2.phi[i],IOB2.replication[i],\
                  tag, recv_rank};
    lc_halo_send2.push_back(h);
    ++lc_halo_send_cnt2[recv_rank];        // count the send items
  }                                       // prepare to send by contiguous destinations

  if (rank == 0)
    cout << "Done packing LC halos" << endl;

  std::sort(lc_halo_send.begin(),lc_halo_send.end(),comp_by_lc_dest);
  std::sort(lc_halo_send2.begin(),lc_halo_send2.end(),comp_by_lc_dest);

  // END LC HALOS

  if (rank == 0)
    cout << "Done sorting LC halos" << endl;

  clear_IO_buffers();
  MPI_Barrier(Partition::getComm());

  MPI_Barrier(Partition::getComm());

  // create the MPI types
  MPI_Datatype lc_halo_type;
  {
    MPI_Datatype type[6] = { MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_INT64_T, MPI_UNSIGNED, MPI_UB };
    int blocklen[6] = {7,1,1,1,1,1};
    MPI_Aint disp[6] = {  offsetof(lc_halo,posvel_a_m),
                          offsetof(lc_halo,phi),
                          offsetof(lc_halo,rr),
                          offsetof(lc_halo,id),
                          offsetof(lc_halo,destination_rank),
                          sizeof(lc_halo) };
    MPI_Type_struct(6,blocklen,disp,type,&lc_halo_type);
    MPI_Type_commit(&lc_halo_type);
  }


  // get the receive counts
  vector<int> lc_halo_recv_cnt;
  lc_halo_recv_cnt.resize(n_ranks,0);
  MPI_Alltoall(&lc_halo_send_cnt[0],1,MPI_INT,&lc_halo_recv_cnt[0],1,MPI_INT,Partition::getComm());

  vector<int> lc_halo_recv_cnt2;
  lc_halo_recv_cnt2.resize(n_ranks,0);
  MPI_Alltoall(&lc_halo_send_cnt2[0],1,MPI_INT,&lc_halo_recv_cnt2[0],1,MPI_INT,Partition::getComm());
  // each rank now knows how many items it will receive from every other rank

  // calculate the offsets
  vector<int> lc_halo_send_off;
  lc_halo_send_off.resize(n_ranks,0);
  vector<int> lc_halo_send_off2;
  lc_halo_send_off2.resize(n_ranks,0);


  vector<int> lc_halo_recv_off;
  lc_halo_recv_off.resize(n_ranks,0);
  vector<int> lc_halo_recv_off2;
  lc_halo_recv_off2.resize(n_ranks,0);

  lc_halo_send_off[0] = lc_halo_recv_off[0] = 0;
  lc_halo_send_off2[0] = lc_halo_recv_off2[0] = 0;

  for (int i=1; i<n_ranks; ++i) {
    lc_halo_send_off[i] = lc_halo_send_off[i-1] + lc_halo_send_cnt[i-1];
    lc_halo_recv_off[i] = lc_halo_recv_off[i-1] + lc_halo_recv_cnt[i-1];
    lc_halo_send_off2[i] = lc_halo_send_off2[i-1] + lc_halo_send_cnt2[i-1];
    lc_halo_recv_off2[i] = lc_halo_recv_off2[i-1] + lc_halo_recv_cnt2[i-1];
  }

  // compute the receive totals to allocate buffers
  int lc_halo_recv_total = 0;
  int lc_halo_recv_total2 = 0;

  for (int i=0; i<n_ranks; ++i) {
    lc_halo_recv_total += lc_halo_recv_cnt[i];
    lc_halo_recv_total2 += lc_halo_recv_cnt2[i];

  }

  vector<lc_halo> lc_halo_recv;
  lc_halo_recv.resize(lc_halo_recv_total);
  vector<lc_halo> lc_halo_recv2;
  lc_halo_recv2.resize(lc_halo_recv_total2);

  MPI_Barrier(Partition::getComm());

  if (rank == 0)
    cout << "About to send data" << endl;


  // send the actual data
  MPI_Alltoallv(&lc_halo_send[0],&lc_halo_send_cnt[0],&lc_halo_send_off[0],lc_halo_type,\
                   &lc_halo_recv[0],&lc_halo_recv_cnt[0],&lc_halo_recv_off[0],lc_halo_type,Partition::getComm());

  MPI_Alltoallv(&lc_halo_send2[0],&lc_halo_send_cnt2[0],&lc_halo_send_off2[0],lc_halo_type,\
                   &lc_halo_recv2[0],&lc_halo_recv_cnt2[0],&lc_halo_recv_off2[0],lc_halo_type,Partition::getComm());

  if (rank == 0)
    cout << "About to sort" << endl;

   std::sort(lc_halo_recv.begin(),lc_halo_recv.end(),comp_by_id);
   std::stable_sort(lc_halo_recv.begin(),lc_halo_recv.end(),comp_by_rep);
   std::sort(lc_halo_recv2.begin(),lc_halo_recv2.end(),comp_by_id);
   std::stable_sort(lc_halo_recv2.begin(),lc_halo_recv2.end(),comp_by_rep);

   // at this point we have saved all the data 
  if (rank == 0)
    cout << "Sorted" << endl;


 clear_IO_buffers();
  for (int i=0; i<lc_halo_recv_total; ++i) {
   // halo_properties_t tmp = fof_lookup(lc_halo_recv[i].id, fof_halo_recv);
    IOB.x.push_back(lc_halo_recv[i].posvel_a_m[0]);
    IOB.y.push_back(lc_halo_recv[i].posvel_a_m[1]);
    IOB.z.push_back(lc_halo_recv[i].posvel_a_m[2]);
    IOB.vx.push_back(lc_halo_recv[i].posvel_a_m[3]);
    IOB.vy.push_back(lc_halo_recv[i].posvel_a_m[4]);
    IOB.vz.push_back(lc_halo_recv[i].posvel_a_m[5]);
    IOB.a.push_back(lc_halo_recv[i].posvel_a_m[6]);
    IOB.replication.push_back(lc_halo_recv[i].rr);
    IOB.id.push_back(lc_halo_recv[i].id);
    IOB.phi.push_back(lc_halo_recv[i].phi);
  }
  for (int i=0; i<lc_halo_recv_total2; ++i) {
   // halo_properties_t tmp = fof_lookup(lc_halo_recv[i].id, fof_halo_recv);
    IOB2.x.push_back(lc_halo_recv2[i].posvel_a_m[0]);
    IOB2.y.push_back(lc_halo_recv2[i].posvel_a_m[1]);
    IOB2.z.push_back(lc_halo_recv2[i].posvel_a_m[2]);
    IOB2.vx.push_back(lc_halo_recv2[i].posvel_a_m[3]);
    IOB2.vy.push_back(lc_halo_recv2[i].posvel_a_m[4]);
    IOB2.vz.push_back(lc_halo_recv2[i].posvel_a_m[5]);
    IOB2.a.push_back(lc_halo_recv2[i].posvel_a_m[6]);
    IOB2.replication.push_back(lc_halo_recv2[i].rr);
    IOB2.id.push_back(lc_halo_recv2[i].id);
    IOB2.phi.push_back(lc_halo_recv2[i].phi);
  }

  if (rank == 0)
    cout << "Assigned to buffers" << endl;


  if (lc_halo_recv_total != lc_halo_recv_total2){
    err += 1;
    //if (rank==0)
      cout << "The number of elements is different in the two files for rank " << rank << endl;
      cout << lc_halo_recv_total << " , and " << lc_halo_recv_total2 << " for rank "<< rank<< endl;
    }
   else{

    for (int i=0;i<lc_halo_recv_total;i++){
       if((IOB.id[i]!=IOB2.id[i])||(IOB.replication[i]!=IOB2.replication[i]))
         err+=1;
     }
   }



    if (skip_err&&(err>0)){
       err = 0;
       cout << "Skipping over missing particles for rank "<< rank << endl;
       int i = 0;
       while (i < IOB.id.size()){
          if ((IOB.id[i]==IOB2.id[i])&&(IOB.replication[i]==IOB2.replication[i])){
             i += 1;
           }
          else {
           bool not_found = true;
           for (int j=0;j<32;j++){
               if ((IOB.id[i]==IOB2.id[i+j])&&(IOB.replication[i]==IOB2.replication[i+j])){
                  for (int k=0; k<j; k++)
                      erase_element_IOB2(i+k);
                     //IOB2.id.erase(IOB2.id.begin()+i+k); // if the particle is found within the next 32 elements then delete section before that
                not_found = false;
                i+=1;
                }
              }
           if (not_found)
             erase_element_IOB(i); 
      //       IOB.id.erase(IOB.id.begin()+i); // delete and just continue until we find 
         }
     }

    cout<< "number of elements before was " <<  lc_halo_recv_total2 << " , and "<<lc_halo_recv_total <<endl;
    cout<< "number of elements now is " <<  IOB.id.size() << " , and "<< IOB2.id.size() <<endl;
    // err = 1; // for now, to supress later outputs
    lc_halo_recv_total = IOB.id.size();
    lc_halo_recv_total2 = IOB2.id.size();
    }


    float da=0;
    float dx=0;
    float dy=0;
    float dz=0;
    float dvx=0;
    float dvy=0;
    float dvz=0;
    int64_t diff_id =0;
    int diff_rep = 0;
    float dphi=0;

//  if (lc_halo_recv_total == lc_halo_recv_total2){
    if (err==0){
    for (int i=0;i<lc_halo_recv_total;++i){
         float xdiff = fabs(IOB.x[i]-IOB2.x[i]);
         float ydiff = fabs(IOB.y[i]-IOB2.y[i]);
         float zdiff = fabs(IOB.z[i]-IOB2.z[i]);
         float vxdiff = fabs(IOB.vx[i]-IOB2.vx[i]);
         float vydiff = fabs(IOB.vy[i]-IOB2.vy[i]);
         float vzdiff = fabs(IOB.vz[i]-IOB2.vz[i]);
         float adiff = fabs(IOB.a[i]-IOB2.a[i]);
         float phidiff = fabs(IOB.phi[i]-IOB2.phi[i]);
         dx = (dx<xdiff)?xdiff:dx;
         dy = (dy<ydiff)?ydiff:dy;
         dz = (dz<zdiff)?zdiff:dz;
         dvx = (dvx<vxdiff)?vxdiff:dvx;
         dvy = (dvy<vydiff)?vydiff:dvy;
         dvz = (dvz<vzdiff)?vzdiff:dvz;
         da = (da<adiff)?adiff:da;
         int64_t  iddiff = abs(IOB.id[i]-IOB2.id[i]);
         int repdiff = abs(IOB.replication[i]-IOB2.replication[i]);
         diff_id = (iddiff<diff_id)?diff_id:iddiff;
         diff_rep = (repdiff<diff_rep)?diff_rep:repdiff;
         dphi = (dphi<phidiff)?phidiff:dphi;
      }
   }

    if ((diff_id>0)||(diff_rep>0)){
      err +=1;
      if (rank==0)
       cout <<"Ids don't match up!"<<endl;
    }




    float dx_tot;
    float dy_tot;
    float dz_tot;
    float dvx_tot;
    float dvy_tot;
    float dvz_tot;
    float da_tot; 
    float dphi_tot; 
    int err_tot;

   if (rank == 0)
    cout << "Starting reduction" << endl;
    MPI_Reduce(&dx, &dx_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&dy, &dy_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&dz, &dz_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&dvx, &dvx_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&dvy, &dvy_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&dvz, &dvz_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&da, &da_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&dphi, &dphi_tot, 1, MPI_FLOAT, MPI_MAX, 0, Partition::getComm());
    MPI_Reduce(&err, &err_tot, 1, MPI_INT, MPI_SUM, 0, Partition::getComm());



  if ((rank == 0) && (err_tot ==0)){
    cout << "Maximum dx = "<< dx_tot << endl;
    cout << "Maximum dy = "<< dy_tot << endl;
    cout << "Maximum dz = "<< dz_tot << endl;
    cout << "Maximum dvx = "<< dvx_tot << endl;
    cout << "Maximum dvy = "<< dvy_tot << endl;
    cout << "Maximum dvz = "<< dvz_tot << endl;
    cout << "Maximum da = "<< da_tot << endl;
    cout << "Maximum dphi = "<< dphi_tot << endl;

  }

// communicate these between ranks and find the maximum

  // now if the difference in values is >0 , then create histograms of the differences to check they're unbiased
  // compute the mean and standard deviation of the differences and look at the fractional errors.   
  if (err_tot == 0){
  if ((dx_tot>0)||(dy_tot>0)||(dz_tot>0)){
    compute_mean_std_dist(IOB.x,IOB2.x,dx_tot,"x position");
    compute_mean_std_dist(IOB.y,IOB2.y,dy_tot,"y position");
    compute_mean_std_dist(IOB.z,IOB2.z,dz_tot,"z position");
  }
  if ((dvx_tot>0)||(dvy_tot>0)||(dvz_tot>0)){
    compute_mean_std_dist(IOB.vx,IOB2.vx,dvx_tot,"x velocity");
    compute_mean_std_dist(IOB.vy,IOB2.vy,dvy_tot,"y velocity");
    compute_mean_std_dist(IOB.vz,IOB2.vz,dvz_tot,"z velocity");
  }
  if (da_tot>0){
    compute_mean_std_dist(IOB.a,IOB2.a,da_tot,"scale factor");
  }
  if (dphi_tot>0){
    compute_mean_std_dist(IOB.phi,IOB2.phi,dphi_tot,"potential");
  }
  }
 


  MPI_Barrier(Partition::getComm());


  Partition::finalize();
  MPI_Finalize();
  return 0;
}
