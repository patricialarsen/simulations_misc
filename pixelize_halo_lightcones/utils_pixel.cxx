#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <getopt.h>

#if defined(_OPENMP) 
#include <omp.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <stdint.h>

#include <iostream>
#include <fstream>

#include <vector>
#include <climits>
#include <algorithm>
#include "chealpix.h"
#include "pointing.h"
#include "healpix_base.h"
#include "vec3.h"



// dtk and pixelize header files
#include "pixelize_halos_sod.h"
#include "hdf5_util.hpp"
#include "GenericIO.h"

#define POSVEL_T float
#define ID_T int64_t
#define MASK_T uint16_t


using namespace std;
using namespace gio;

MPI_Datatype createparticles()
{
    MPI_Datatype particles_mpi;
    MPI_Datatype type[17] = {MPI_FLOAT,MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,MPI_INT,MPI_INT, MPI_INTEGER8 };
    int blocklen[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    MPI_Aint disp[17] = {offsetof(particles,x),
                        offsetof(particles,y),
                        offsetof(particles,z),
                        offsetof(particles,vx),
                        offsetof(particles,vy),
                        offsetof(particles,vz),
                        offsetof(particles,rot),
                        offsetof(particles,rep),
                        offsetof(particles,fof_mass),
                        offsetof(particles,a),
                        offsetof(particles,sod_mass),
                        offsetof(particles,sod_radius),
                        offsetof(particles,sod_cdelta),
                        offsetof(particles,sod_cdelta_error),
                        offsetof(particles,rank),
                        offsetof(particles,pix_index),
                        offsetof(particles,id)};
    MPI_Type_struct(17,blocklen,disp,type,&particles_mpi);
    MPI_Type_commit(&particles_mpi);
    return particles_mpi;
}


bool comp_rank(const particles &a, const particles &b){
   return a.rank < b.rank;
}


int compute_rank_fullsky(int pix_val, int numranks){
    int rankn = pix_val%numranks;
}

int compute_rank_partsky(int pix_val, int numranks, vector<int> pix_list){
   // assign to MPI ranks (partial-sky version)
   int rankn=-99;
   int i;
   for (i=0;i<pix_list.size();++i){
       if (pix_list[i]==pix_val){
       rankn = i%numranks;
       return rankn;
      }
   }
  if (rankn==-99){
   printf("pix_val = %d \n", pix_val);
   printf("WARNING - OCTANT CHOSEN WHEN PARTICLES OUTSIDE OF THIS AREA ARE PRESENT, COMPUTING FULLSKY VERSION \n");
   }
   rankn = compute_rank_fullsky(pix_val,numranks);
   return rankn;
}

int compute_oct_pixels(T_Healpix_Base<int> map_hires , T_Healpix_Base<int> map_lores, vector<int> &pix_list_oct, vector<int> &pix_list_noct, double edge_val){
    long npix1 = map_lores.Npix();
    long npix = map_hires.Npix();

    for (int k=0;k<npix;++k){
        vec3 pix_center = map_hires.pix2vec(k);
        if ((pix_center.x >= -edge_val) && (pix_center.y>=-edge_val) && (pix_center.z<=edge_val)){
            int ind_lo = map_lores.vec2pix(pix_center);
            int check_pix = 0;
            for (int i=0;i<pix_list_oct.size();++i){
                if (ind_lo==pix_list_oct[i]){
                    check_pix +=1;
                }
            }
            if (check_pix==0){
                pix_list_oct.push_back(ind_lo);
            }
        }
    }
    for (int i=0;i<npix1;++i){
        int check_pix = 0;
        for (int k=0;k<pix_list_oct.size();++k){
            if (pix_list_oct[k]==i){
                check_pix +=1;
            }
        }
        if (check_pix==0){
            pix_list_noct.push_back(i);
        }
    }
    return 0;
}


int read_file(vector<POSVEL_T> &xx, vector<POSVEL_T> &yy, vector<POSVEL_T> &zz, vector<POSVEL_T> &vx, vector<POSVEL_T> &vy, vector<POSVEL_T> &vz, vector<POSVEL_T> &a, vector<ID_T> &id, vector<int> &rotation, 
    vector<int> &replication, vector<POSVEL_T> &fof_mass, vector<POSVEL_T> &sod_mass, vector<POSVEL_T> &sod_radius, vector<POSVEL_T> &sod_cdelta, vector<POSVEL_T> &sod_cdelta_error, size_t &Np, string mpiioName){
    unsigned Method = GenericIO::FileIOPOSIX;
    const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
    printf("filename_read = %s \n", mpiioName.c_str());
    if (EnvStr && string(EnvStr) == "1")
        Method = GenericIO::FileIOMPI;
        GenericIO GIO( MPI_COMM_WORLD, mpiioName, Method);
        {
        GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
        MPI_Barrier(MPI_COMM_WORLD);
        
        Np = GIO.readNumElems();
        printf("Np = %lu \n", Np);

        xx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        yy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        zz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        vx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        vy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        vz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        rotation.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
        replication.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
        fof_mass.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        sod_mass.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        sod_cdelta.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        sod_cdelta_error.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        sod_radius.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));

        a.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
        id.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));

        GIO.addVariable("x", xx, true);
        GIO.addVariable("y", yy, true);
        GIO.addVariable("z", zz, true);
        GIO.addVariable("vx", vx, true);
        GIO.addVariable("vy", vy, true);
        GIO.addVariable("vz", vz, true);

        GIO.addVariable("rotation", rotation, true);
        GIO.addVariable("replication", replication, true);
        GIO.addVariable("fof_mass", fof_mass, true);
        GIO.addVariable("sod_mass", sod_mass, true);
        GIO.addVariable("sod_radius", sod_radius, true);
        GIO.addVariable("sod_cdelta", sod_cdelta, true);
        GIO.addVariable("sod_cdelta_error", sod_cdelta_error, true);

        GIO.addVariable("a", a, true);
        GIO.addVariable("id", id, true);

        GIO.readData();
        }
    return 0;
}

int redistribute_particles(int numranks,vector<particles> &particle_data1, vector<particles> &recv_particles, MPI_Datatype particles_mpi, vector<int> &send_count, int &nparticles){
    vector<int> send_offset;
    send_offset.resize(numranks);
    send_offset[0] = 0;

    vector<int> recv_offset;
    recv_offset.resize(numranks);
    recv_offset[0] = 0;

    vector<int> recv_count;
    recv_count.resize(numranks);


    MPI_Alltoall( &send_count[0], 1 , MPI_INT, &recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);

    for (int k=1;k<numranks;k++){
        send_offset[k] = send_offset[k-1] + send_count[k-1];
        recv_offset[k] = recv_offset[k-1] + recv_count[k-1];
    }
    recv_particles.resize(recv_offset.back()+recv_count.back());

    MPI_Alltoallv(&particle_data1[0], &send_count[0], &send_offset[0],particles_mpi,\
                  &recv_particles[0], &recv_count[0], &recv_offset[0], particles_mpi, MPI_COMM_WORLD );

    nparticles = recv_offset.back()+recv_count.back();
    printf("nparticles = %d\n",nparticles);

    return 0;
}




int compute_ranks(int Np, vector<POSVEL_T> &xx, vector<POSVEL_T> &yy, vector<POSVEL_T> &zz, vector<POSVEL_T> &vx, vector<POSVEL_T> &vy, vector<POSVEL_T> &vz, vector<POSVEL_T> &a, vector<ID_T> &id, vector<int> &rot,
     vector<int> &rep, vector<POSVEL_T> &fof_mass, vector<POSVEL_T> &sod_mass, vector<POSVEL_T> &sod_radius, vector<POSVEL_T> &sod_cdelta, vector<POSVEL_T> &sod_cdelta_error, T_Healpix_Base<int> map_lores, 
     vector<int> pix_list_oct, int numranks, vector<particles> &particle_data1, vector<int> &pix_id, int octant, vector<int> &send_count, string step_s){
   // compute ranks for halos
   for (int k=0;k<Np;k++){
     int dist_check2 = xx[k]*xx[k] + yy[k]*yy[k] + zz[k]*zz[k];
     if (dist_check2<25000000){
       zz[k] = -zz[k];
       int rankn;
       vec3 vec_temp3 = vec3(xx[k],yy[k],zz[k]);
       int ind2 = map_lores.vec2pix(vec_temp3);
       pix_id[k] = ind2;
       rankn = compute_rank_partsky(ind2, numranks, pix_list_oct);
       particles part = {xx[k], yy[k], zz[k], vx[k],vy[k],vz[k], rot[k], rep[k], fof_mass[k], a[k], sod_mass[k], sod_radius[k], sod_cdelta[k], sod_cdelta_error[k], rankn, ind2 , id[k]};
       if (step_s=="347"){
         float a_cut = 0.69751245;
         if (a[k]>a_cut){
           particle_data1.push_back(part);
           send_count[rankn] += 1;
         }
       }
       else if (step_s=="338"){
         float a_cut = 0.69751245;
           if (a[k]<=a_cut){
             particle_data1.push_back(part);
             send_count[rankn] += 1;
           }
         }
         else if (step_s=="259"){
           float a_cut = 0.52238804;
           if (a[k]>a_cut){
           send_count[rankn] += 1;
           particle_data1.push_back(part);
         }
       }
       else if (step_s=="253"){
         float a_cut = 0.52238804;
         if (a[k]<=a_cut){
           send_count[rankn] += 1;
           particle_data1.push_back(part);
         }
       }
       else{
         send_count[rankn] += 1;
         particle_data1.push_back(part);
       }
     }
   }
 return 0;
}



int read_and_redistribute(string name, T_Healpix_Base<int> map_lores, int numranks, int octant, vector<int> &send_count, vector<particles> &recv_particles, MPI_Datatype particles_mpi, int &nparticles, vector<int> &pix_list_oct, string step_s){
    vector<particles> particle_data1;
    {
        vector<POSVEL_T> xx, yy, zz;
        vector<POSVEL_T> vx, vy, vz;
        vector<int> rot, rep;
        vector<POSVEL_T> fof_mass;
        vector<POSVEL_T> sod_mass, sod_radius, sod_cdelta, sod_cdelta_error;
        vector<POSVEL_T> a;
        vector<int> pix_id;
        assert(sizeof(ID_T) == 8);
        vector<ID_T> id;

        size_t Np = 0;
        int status = read_file(xx,yy,zz,vx,vy,vz,a,id,rot,rep,fof_mass,sod_mass,sod_radius,sod_cdelta,sod_cdelta_error,Np,name);
        if(Np!=0){
            printf("Np>0\n");
            xx.resize(Np);
            yy.resize(Np);
            zz.resize(Np);
            a.resize(Np);
            pix_id.resize(Np);
            id.resize(Np);
            fof_mass.resize(Np);
            sod_mass.resize(Np);

            sod_radius.resize(Np);
            sod_cdelta.resize(Np);
            sod_cdelta_error.resize(Np);

            vx.resize(Np);
            vy.resize(Np);
            vz.resize(Np);
            rot.resize(Np);
            rep.resize(Np);

           status =  compute_ranks(Np, xx, yy, zz, vx, vy, vz, a, id, rot, rep, fof_mass, sod_mass, sod_radius, sod_cdelta, sod_cdelta_error, map_lores,
               pix_list_oct, numranks, particle_data1, pix_id,  octant, send_count, step_s);
 
           sort(particle_data1.begin(),particle_data1.end(),comp_rank);

        }
    }

    int status = redistribute_particles(numranks, particle_data1,  recv_particles,   particles_mpi, send_count, nparticles);
    return status;
}



int write_files(vector<particles> recv_particles, int nparticles, int l, vector<int> pixel_ids, string zrange_str, string step_s){

        stringstream ss;
        ss<< pixel_ids[l];
        string str_pix = ss.str();

        H5::H5File file("/projects/DarkUniverse_esp/prlarsen/healpix_cutouts_dc2_32_5000run_test/" + zrange_str + "/cutout_" + str_pix + ".hdf5", H5F_ACC_RDWR);
        string file_name1 = "/projects/DarkUniverse_esp/prlarsen/healpix_cutouts_dc2_32_5000run_test/" + zrange_str + "/cutout_" + str_pix + ".hdf5";


   vector<POSVEL_T> x_1;
   vector<POSVEL_T> y_1;
   vector<POSVEL_T> z_1;
   vector<POSVEL_T> vx_1;
   vector<POSVEL_T> vy_1;
   vector<POSVEL_T> vz_1;

   vector<int> rot_1;
   vector<int> rep_1;
   vector<POSVEL_T> a_1;
   vector<ID_T>  id_1;
   vector<POSVEL_T> fof_mass_1, sod_mass_1;
   vector<POSVEL_T> sod_radius_1, sod_cdelta_1, sod_cdelta_error_1;

   for (int i=0;i<nparticles;++i){

   if (recv_particles[i].pix_index==pixel_ids[l]){
     x_1.push_back(recv_particles[i].x);
     y_1.push_back(recv_particles[i].y);
     z_1.push_back(recv_particles[i].z);
     vx_1.push_back(recv_particles[i].vx);
     vy_1.push_back(recv_particles[i].vy);
     vz_1.push_back(recv_particles[i].vz);
     a_1.push_back(recv_particles[i].a);
     fof_mass_1.push_back(recv_particles[i].fof_mass);
     sod_mass_1.push_back(recv_particles[i].sod_mass);
     sod_radius_1.push_back(recv_particles[i].sod_radius);
     sod_cdelta_1.push_back(recv_particles[i].sod_cdelta);
     sod_cdelta_error_1.push_back(recv_particles[i].sod_cdelta_error);
     rot_1.push_back(recv_particles[i].rot);
     rep_1.push_back(recv_particles[i].rep);
     id_1.push_back(recv_particles[i].id);
   }
   }

   // variable names 
   hsize_t x1size = x_1.size();
   std::string var_name_x, var_name_y, var_name_z;
   std::string var_name_vx, var_name_vy, var_name_vz;
   std::string var_name_id, var_name_a, var_name_rot, var_name_rep;
   std::string var_name_fof_mass, var_name_sod_mass, var_name_sod_radius;
   std::string var_name_sod_cdelta, var_name_sod_cdelta_error;


   var_name_x="/" + step_s + "/x";
   var_name_y="/" + step_s + "/y";
   var_name_z="/" + step_s + "/z";
   var_name_vx="/" + step_s + "/vx";
   var_name_vy="/" + step_s + "/vy";
   var_name_vz="/" + step_s + "/vz";
   var_name_rep="/" + step_s + "/rep";
   var_name_rot="/" + step_s + "/rot";
   var_name_a="/" + step_s + "/a";
   var_name_id="/" + step_s + "/id";
   var_name_fof_mass="/"+step_s+"/fof_mass";
   var_name_sod_mass="/"+step_s+"/sod_mass";
   var_name_sod_radius="/"+step_s+"/sod_radius";
   var_name_sod_cdelta="/"+step_s+"/sod_cdelta";
   var_name_sod_cdelta_error="/"+step_s+"/sod_cdelta_error";



   dtk::write_hdf5(file, var_name_x,  &x_1[0] ,x1size);
   dtk::write_hdf5(file, var_name_y,  &y_1[0] ,x1size);
   dtk::write_hdf5(file, var_name_z,  &z_1[0] ,x1size);
   dtk::write_hdf5(file, var_name_vx,  &vx_1[0] ,x1size);
   dtk::write_hdf5(file, var_name_vy,  &vy_1[0] ,x1size);
   dtk::write_hdf5(file, var_name_vz,  &vz_1[0] ,x1size);
   dtk::write_hdf5(file, var_name_rot, &rot_1[0], x1size);
   dtk::write_hdf5(file, var_name_rep, &rep_1[0], x1size);
   dtk::write_hdf5(file, var_name_a, &a_1[0], x1size);
   dtk::write_hdf5(file, var_name_id, &id_1[0], x1size);
   dtk::write_hdf5(file, var_name_fof_mass, &fof_mass_1[0], x1size);
   dtk::write_hdf5(file, var_name_sod_mass, &sod_mass_1[0], x1size);
   dtk::write_hdf5(file, var_name_sod_radius, &sod_radius_1[0], x1size);
   dtk::write_hdf5(file, var_name_sod_cdelta, &sod_cdelta_1[0], x1size);
   dtk::write_hdf5(file, var_name_sod_cdelta_error, &sod_cdelta_error_1[0], x1size);

  return 0;
}
