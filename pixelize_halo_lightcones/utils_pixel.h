#define POSVEL_T float
#define ID_T int64_t
#include <vector>

using namespace std;

MPI_Datatype createparticles();

int compute_rank_fullsky(int pix_val, int numranks);

int compute_rank_partsky(int pix_val, int numranks, vector<int> pix_list);

int compute_oct_pixels(T_Healpix_Base<int> map_hires , T_Healpix_Base<int> map_lores, vector<int> &pix_list_oct, vector<int> &pix_list_noct, double edge_val);

int read_file(vector<POSVEL_T> &xx, vector<POSVEL_T> &yy, vector<POSVEL_T> &zz, vector<POSVEL_T> &vx, vector<POSVEL_T> &vy, vector<POSVEL_T> &vz, vector<POSVEL_T> &a, vector<ID_T> &id, vector<int> &rotation, vector<int> &replication, 
    vector<POSVEL_T> &fof_mass, vector<POSVEL_T> &sod_mass, vector<POSVEL_T> &sod_radius, vector<POSVEL_T> &sod_cdelta, vector<POSVEL_T> &sod_cdelta_error, size_t &Np, string mpiioName);

int compute_ranks(int Np, vector<POSVEL_T> &xx, vector<POSVEL_T> &yy, vector<POSVEL_T> &zz, vector<POSVEL_T> &vx, vector<POSVEL_T> &vy, vector<POSVEL_T> &vz, vector<POSVEL_T> &a, vector<ID_T> &id, vector<int> &rot,
     vector<int> &rep, vector<POSVEL_T> &fof_mass, vector<POSVEL_T> &sod_mass, vector<POSVEL_T> &sod_radius, vector<POSVEL_T> &sod_cdelta, vector<POSVEL_T> &sod_cdelta_error, T_Healpix_Base<int> map_lores,
     vector<int> pix_list_oct, int numranks, vector<particles> &particle_data1, vector<int> &pix_id, int octant, vector<int> &send_count, string step_s);


int read_and_redistribute(string name, T_Healpix_Base<int> map_lores, int numranks, int octant, vector<int> &send_count, vector<particles> &recv_particles, MPI_Datatype particles_mpi, int &nparticles, vector<int> &pix_list_oct, string step_s);

int write_files(vector<particles> recv_particles, int nparticles, int l, vector<int> pixel_ids, string zrange_str, string step_s);

