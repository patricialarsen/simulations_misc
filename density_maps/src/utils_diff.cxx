26d25
< #include <fstream>
29d27
< #include <random>
43,46c41,42
< #include <filesystem>
< 
< #include "utils.h"
< #include "pix_funcs.h"
---
> #include "utils_gravonly.h"
> #include "pix_funcs_gravonly.h"
50a47
> 
58a56
> 
63d60
< 
78,79d74
< 
< 
82c77
<   GenericIO GIO(MPI_COMM_WORLD,file_name);
---
>   GenericIO GIO(MPI_COMM_WORLD,file_name,GenericIO::FileIOMPI);
101a97
> 
108c104
< 	  
---
> 
150c146
<   //Map to rank based on id. Looks complicated im just jumbling id (using Thomas Wang's 64-bit integer hash function) so its balanced. 
---
>     //Map to rank based on id. Looks complicated im just jumbling id (using Thomas Wang's 64-bit integer hash function) so its balanced.
172d167
< 
206d200
< 
242d235
< 
309c302
<   double PhysOrigin[3]; 
---
>   double PhysOrigin[3];
325,326c318,319
< 	  NewGIO.setPhysOrigin(PhysOrigin[d], d);
< 	  NewGIO.setPhysScale(PhysScale[d], d);
---
>           NewGIO.setPhysOrigin(PhysOrigin[d], d);
>           NewGIO.setPhysScale(PhysScale[d], d);
333,340d325
< 
<   //NOTE: warning here that random_shuffle is deprecated because of default behaviour. 
<   // I think it should be fine but support is limited
<   random_device rd; 
<   mt19937 g(rd());
<   shuffle(indices.begin(), indices.end(),g);
<   //random_shuffle(indices.begin(), indices.end(), drand48elmt);
< 
466c451,452
< void read_and_redistribute(string file_name, string file_name_next, int numranks, PLParticles* P,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff, bool output_downsampled, float downsampling_rate, string file_name_output){
---
> void read_and_redistribute(string file_name, string file_name_next, int numranks, PLParticles* P,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff,bool output_downsampled, float downsampling_rate, string file_name_output){
> 
471a458
> 	
478c465
<     read_particles(P, file_name, file_name_next);    
---
>     read_particles(P, file_name, file_name_next);
490,491c477,478
< 	assert(downsampling_rate>0.0);
< 	output_downsampled_particles( P, downsampling_rate, file_name_output, file_name);
---
>         assert(downsampling_rate>0.0);
>         output_downsampled_particles( P, downsampling_rate, file_name_output, file_name);
502c489
<     status = compute_ranks_count( P,  map_lores, map_hires, numranks, send_count, rank_diff); 
---
>     status = compute_ranks_count( P,  map_lores, map_hires, numranks, send_count, rank_diff);
519d505
< 
539,542d524
<     if (start_idx.size()==0){
<       MPI_File_iwrite_at_all(fh,0,NULL,0,MPI_DOUBLE, &req);
<       MPI_Wait(&req, MPI_STATUS_IGNORE);
<     }	    
547c529,530
<     MPI_File_iwrite_at_all(fh,offset,&rho[start_tmp],off_tmp,MPI_DOUBLE, &req);
---
>     MPI_File_seek(fh,offset,MPI_SEEK_SET);
>     MPI_File_iwrite(fh,&rho[start_tmp],off_tmp,MPI_DOUBLE, &req);
553d535
< void write_files_hydro(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<double> rho, vector<double> phi, vector<double> vel, vector<double> ksz, vector<double> tsz, vector<double> xray1, vector<double> xray2, vector<double> xray3, vector<double> xray4, vector<double> temp, int64_t npix_hires){
554a537,538
> 
> void write_files(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<double> rho, vector<double> phi, vector<double> vel, int64_t npix_hires){
594,596c578,579
<   rho.resize(size_pad); phi.resize(size_pad); vel.resize(size_pad); ksz.resize(size_pad); tsz.resize(size_pad); 
<   xray1.resize(size_pad); xray2.resize(size_pad); xray3.resize(size_pad); xray4.resize(size_pad);
<   temp.resize(size_pad); data_idx.resize(size_pad);
---
>   rho.resize(size_pad); phi.resize(size_pad); vel.resize(size_pad); 
>   data_idx.resize(size_pad);
601,607d583
<   GIO.addVariable("ksz", ksz, GenericIO::VarHasExtraSpace);
<   GIO.addVariable("tsz", tsz, GenericIO::VarHasExtraSpace);
<   GIO.addVariable("xray1", xray1, GenericIO::VarHasExtraSpace);//ROSAT
<   GIO.addVariable("xray2", xray2, GenericIO::VarHasExtraSpace);//ErositaLo
<   GIO.addVariable("xray3", xray3, GenericIO::VarHasExtraSpace);//ErositaHi
<   GIO.addVariable("xray4", xray4, GenericIO::VarHasExtraSpace);//Bolo
<   GIO.addVariable("temp", temp, GenericIO::VarHasExtraSpace);
613,655c589,591
<   rho.resize(size_n); phi.resize(size_n); vel.resize(size_n); ksz.resize(size_n); tsz.resize(size_n); 
<   xray1.resize(size_n); xray2.resize(size_n); xray3.resize(size_n); xray4.resize(size_n);
<   temp.resize(size_n); data_idx.resize(size_n);
<   }
< 
< 
< 
< void write_files(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<float> rho, vector<float> phi, vector<float> vel, int64_t npix_hires){ 
< 
<   int commRank, commRanks;
<   MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
<   MPI_Comm_size(MPI_COMM_WORLD, &commRanks);
< 
<   string output_name = outfile + stepnumber + "_dens.bin";
<   string output_name_phi = outfile + stepnumber + "_phi.bin";
<   string output_name_vel = outfile + stepnumber + "_vel.bin";
< 
< 
<   const char *name_out = output_name.c_str();
<   const char *name_out_phi = output_name_phi.c_str();
<   const char *name_out_vel = output_name_vel.c_str();
< 
<   MPI_File fh, fh_phi, fh_vel;
<   MPI_Request req, req_phi, req_vel;
< 
<   MPI_File_open(MPI_COMM_WORLD, name_out, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
<   MPI_File_open(MPI_COMM_WORLD, name_out_phi, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_phi);
<   MPI_File_open(MPI_COMM_WORLD, name_out_vel, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_vel);
< 
<   MPI_Offset filesize = npix_hires*sizeof(float);
<   MPI_File_set_size(fh, filesize);
<   MPI_File_set_size(fh_phi, filesize);
<   MPI_File_set_size(fh_vel, filesize); 
< 
< 
<   int status;
<   status = output_file(commRank, fh, req, rho, start_idx, end_idx,  pix_nums_start);
<   status = output_file(commRank, fh_phi, req_phi, phi, start_idx, end_idx, pix_nums_start);
<   status = output_file(commRank, fh_vel, req_vel, vel, start_idx, end_idx, pix_nums_start);
< 
<   MPI_File_close(&fh);
<   MPI_File_close(&fh_phi);
<   MPI_File_close(&fh_vel);
---
>   rho.resize(size_n); phi.resize(size_n); vel.resize(size_n); 
>   data_idx.resize(size_n);
> }
657,658d592
<   MPI_Barrier(MPI_COMM_WORLD);
<   }
