#include "healpix_base.h"
#include <unordered_map>
#include <string>
#include <vector>
using namespace std;

void initialize_pixel(int pix_val, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, std::vector<float> &rho , std::vector<float> &phi , std::vector<float> &vel, int64_t &count, std::vector<int64_t> &start_idx, std::vector<int64_t> &end_idx, std::vector<int64_t> &pixnum_start, std::vector<int64_t> &pixnum_end, int rank_diff,  unordered_map<int64_t, int64_t>* ring_to_idx);
void get_pix_list_rank(int octant, int rank, int numranks, int64_t npix_lores, std::vector<int> pix_list_oct, std::vector<int> pix_list_noct,  std::vector<int> &lores_pix);
int compute_ranks_count( Particles* P, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int numranks, vector<int> &send_count, int rank_diff);
void compute_ranks_index( Particles* P, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int numranks, vector<int> send_off, vector<int64_t> &id_particles, int rank_diff);
int assign_dm_ngp(vector<float> &rho, vector<float> &phi, vector<float> &vel, Particles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx);
int assign_dm_cic(vector<float> &rho, vector<float> &phi, vector<float> &vel, Particles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx);

int assign_sz_ngp(vector<float> &rho, vector<float> &phi, vector<float> &vel,vector<float> &ksz, vector<float> &tsz, Particles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx);

int assign_sz_cic(vector<float> &rho, vector<float> &phi, vector<float> &vel,vector<float> &ksz, vector<float> &tsz, Particles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx);

void initialize_pixel_hydro(int pix_val,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, vector<float> &rho , vector<float> &phi, vector<float> &vel, vector<float> &ksz, vector<float> &tsz, int64_t &count, vector<int64_t> &start_idx, vector<int64_t> &end_idx, vector<int64_t> &pixnum_start, vector<int64_t> &pixnum_end, int rank_diff,  unordered_map<int64_t, int64_t>* ring_to_idx);
