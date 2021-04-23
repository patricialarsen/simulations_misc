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

#include "LC_test.h"

#include "MurmurHashNeutral2.cpp" 


/// Assumes ids and replication values are consistent but ordering may have changed
/// Redistributes amongst ranks and sorts the values for a one-to-one check of the changes in 
/// each of the main lightcone outputs 


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

LC_test L_1;
LC_test L_2;

/**
typedef struct lc_halo {
  float posvel_a_m[7];
  float phi;
  int rr;
  int64_t id;
  unsigned int destination_rank;
} lc_halo;
*/

bool comp_by_lc_dest(const lc_properties_test &a, const lc_properties_test &b){
  return a.rank < b.rank;
}

/*
bool comp_by_lc_dest(const lc_halo &a, const lc_halo &b) {
  return a.destination_rank < b.destination_rank;
}



bool comp_by_id(const lc_halo &a, const lc_halo &b) {
  return a.id < b.id;
}

bool comp_by_rep(const lc_halo &a, const lc_halo &b) {
  return a.rr < b.rr;
}

*/

bool comp_by_id(const lc_properties_test &a, const lc_properties_test &b){
  return a.id < b.id;
}
bool comp_by_rep(const lc_properties_test &a, const lc_properties_test &b){
  return a.replication < b.replication;
}


int  compute_mean_std_dist(vector<float> *val1 , vector<float> *val2 , int count ,string var_name, float lim){
  // compute the mean and std of the differences and make a histogram to look more closely
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  double diff=0; 
  double diff_frac=0;
  int n_tot;
  double frac_max;
  double mean;
  double stddev=0;

  for (int i=0; i<count; i++){
      diff += (double)(val1->at(i)-val2->at(i));
      if (val1->at(i)!=0){
        double frac = (double)(fabs(val1->at(i)-val2->at(i))/fabs(val1->at(i)));
         diff_frac = (diff_frac<frac)?frac:diff_frac;
      }
   }

   MPI_Allreduce(&diff, &mean, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   MPI_Reduce(&diff_frac, &frac_max, 1, MPI_DOUBLE, MPI_MAX, 0, Partition::getComm());
   MPI_Allreduce(&count, &n_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());  
   mean = mean/n_tot;   

   for (int i=0; i< count; i++){
      stddev += (double) ((val1->at(i)-val2->at(i)-(float)mean)*(val1->at(i)-val2->at(i)-(float)mean)/(n_tot-1));
   }
   double stddev_tot;
   MPI_Reduce(&stddev, &stddev_tot, 1, MPI_DOUBLE, MPI_SUM, 0, Partition::getComm());
   stddev_tot = sqrt(stddev_tot);

   if (rank==0){
     cout << " " << var_name << endl;
     cout << " ______________________________________" <<endl;
     cout << " mean difference = "<< mean << endl;
     cout << " maximum fractional difference = "<< frac_max<< endl;
     cout << " standard deviation of difference = " << stddev_tot << endl;
     cout << endl;
   }
  int err = 0;


  return err;
}


int  compute_mean(LC_test L_1 , LC_test L_2, float lim ){
  // compute the mean and std of the differences and make a histogram to look more closely
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = L_1.num_parts;

  int err=0;
  for (int i =0; i<N_LC_FLOATS; i++){
    string var_name = float_var_names_test[i];
    err += compute_mean_std_dist(L_1.float_data[i],L_2.float_data[i],count,var_name,lim);
  }
  return err;
}



/*
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
*/
inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks) {
    return MurmurHashNeutral2((void*)(&fof_tag),sizeof(int64_t),0) % n_ranks;
}


void read_lc_file(LC_test &L0, string file_name) {
  // reads in relevant lightcone information from a LC particle file
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();
  
  L0.Resize(num_elems + GIO.requestedExtraSpace());
 

  GIO.addVariable("id", *(L0.id), true);
  GIO.addVariable("replication", *(L0.replication), true);
  
  for (int i=0; i<N_LC_FLOATS; ++i)
	  GIO.addVariable((const string) float_var_names_test[i], *(L0.float_data[i]),true);

  GIO.readData();

  L0.Resize(num_elems);
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
  stringstream thresh{ argv[3] };
  float lim{};
  if (!(thresh >> lim))
	  lim = 0.01;


  L_1.Allocate();
  L_2.Allocate();
  L_1.Set_MPIType();
  L_2.Set_MPIType();



  int err = 0;
  bool skip_err = true;

  // LC HALOS
  read_lc_file(L_1, lc_file);                  // read the file
  read_lc_file(L_2, lc_file2);                  // read the file

  if (rank == 0)
    cout << "Done reading LC" << endl;

  vector<lc_properties_test> lc_halo_send;
  vector<int> lc_halo_send_cnt;
  lc_halo_send_cnt.resize(n_ranks,0);
  for (size_t i=0;i<L_1.num_parts;++i) {     // pack a vector of structures
	lc_properties_test tmp = L_1.GetProperties(i);
	tmp.rank = tag_to_rank(tmp.id,n_ranks);
	lc_halo_send.push_back(tmp);
	++lc_halo_send_cnt[tmp.rank];
  }                                       // prepare to send by contiguous destinations
  vector<lc_properties_test> lc_halo_send2;
  vector<int> lc_halo_send_cnt2;
  lc_halo_send_cnt2.resize(n_ranks,0);
  for (size_t i=0;i<L_2.num_parts;++i) {     // pack a vector of structures
        lc_properties_test tmp = L_2.GetProperties(i);
        tmp.rank = tag_to_rank(tmp.id,n_ranks);
        lc_halo_send2.push_back(tmp);
        ++lc_halo_send_cnt2[tmp.rank];
  }                                       // prepare to send by contiguous destinations


  if (rank == 0)
    cout << "Done packing LC buffers" << endl;

  std::sort(lc_halo_send.begin(),lc_halo_send.end(),comp_by_lc_dest);
  std::sort(lc_halo_send2.begin(),lc_halo_send2.end(),comp_by_lc_dest);

  // END LC HALOS

  if (rank == 0)
    cout << "Done sorting LC halos" << endl;

  L_1.Resize(0);
  L_2.Resize(0);



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

  vector<lc_properties_test> lc_halo_recv;
  lc_halo_recv.resize(lc_halo_recv_total);
  vector<lc_properties_test> lc_halo_recv2;
  lc_halo_recv2.resize(lc_halo_recv_total2);

  MPI_Barrier(Partition::getComm());

  if (rank == 0)
    cout << "About to send data" << endl;


  // send the actual data
  MPI_Alltoallv(&lc_halo_send[0],&lc_halo_send_cnt[0],&lc_halo_send_off[0],L_1.LC_properties_MPI_Type,\
                   &lc_halo_recv[0],&lc_halo_recv_cnt[0],&lc_halo_recv_off[0],L_1.LC_properties_MPI_Type,Partition::getComm());

  MPI_Alltoallv(&lc_halo_send2[0],&lc_halo_send_cnt2[0],&lc_halo_send_off2[0],L_2.LC_properties_MPI_Type,\
                   &lc_halo_recv2[0],&lc_halo_recv_cnt2[0],&lc_halo_recv_off2[0],L_2.LC_properties_MPI_Type,Partition::getComm());

  if (rank == 0)
    cout << "About to sort" << endl;

   std::sort(lc_halo_recv.begin(),lc_halo_recv.end(),comp_by_id);
   std::stable_sort(lc_halo_recv.begin(),lc_halo_recv.end(),comp_by_rep);
   std::sort(lc_halo_recv2.begin(),lc_halo_recv2.end(),comp_by_id);
   std::stable_sort(lc_halo_recv2.begin(),lc_halo_recv2.end(),comp_by_rep);

   // at this point we have saved all the data 
  if (rank == 0)
    cout << "Sorted" << endl;

  L_1.Resize(0);
  L_2.Resize(0);

  for (int i=0; i<lc_halo_recv_total; ++i){
     lc_properties_test tmp = lc_halo_recv[i];
     L_1.PushBack(tmp);
  }
    for (int i=0; i<lc_halo_recv_total2; ++i){
     lc_properties_test tmp = lc_halo_recv2[i];
     L_2.PushBack(tmp);
  }
    lc_halo_recv.resize(0);
    lc_halo_recv2.resize(0);


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
       if( (L_1.id->at(i)!=L_2.id->at(i) )|| (L_1.replication->at(i)!=L_2.replication->at(i) ) )
         err+=1;
     }
   }

    int numl1 = L_1.num_parts;
    int numl2 = L_2.num_parts;
    int dn_1 = 0;
    int dn_2 = 0;


    if (skip_err&&(err>0)){
       err = 0;
       cout << "Skipping over missing particles for rank "<< rank << endl;
       int i = 0;
       while (i < numl1){
          if ((L_1.id->at(i)==L_2.id->at(i))&&(L_1.replication->at(i)==L_2.replication->at(i))){
             i += 1;
           }
          else {
           bool not_found = true;
           for (int j=0;j<32;j++){
               if ((L_1.id->at(i)==L_2.id->at(i+j))&&(L_1.replication->at(i)==L_2.replication->at(i+j))){
                  for (int k=0; k<j; k++)
                      L_2.Erase(i+k);
                     //IOB2.id.erase(IOB2.id.begin()+i+k); // if the particle is found within the next 32 elements then delete section before that
                not_found = false;
                i+=1;
                }
              }
           if (not_found){
             L_1.Erase(i); 
	     -- numl1;
	   }
      //       IOB.id.erase(IOB.id.begin()+i); // delete and just continue until we find 
         }
     }

    //cout<< "number of elements before was " <<  lc_halo_recv_total2 << " , and "<<lc_halo_recv_total <<endl;
    //cout<< "number of elements now is " <<  IOB.id.size() << " , and "<< IOB2.id.size() <<endl;
  //  // err = 1; // for now, to supress later outputs
   // lc_halo_recv_total = IOB.id.size();
   // lc_halo_recv
   //
   //_total2 = IOB2.id.size();
    }

    err = compute_mean(L_1,L_2,lim);

/*
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
*/
//  if (lc_halo_recv_total == lc_halo_recv_total2){
  
  
  /*  if (err==0){
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
*/
// communicate these between ranks and find the maximum

  // now if the difference in values is >0 , then create histograms of the differences to check they're unbiased
  // compute the mean and standard deviation of the differences and look at the fractional errors.   
 
 /* if (err_tot == 0){
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
  */
 
 L_1.Deallocate();
 L_2.Deallocate();

  MPI_Barrier(Partition::getComm());


  Partition::finalize();
  MPI_Finalize();
  return 0;
}
