#include "healpix_base.h"
#include <unordered_map>
#include <string>
#include <vector>
using namespace std;

void initialize_pixel(int pix_val, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, std::vector<double> &rho , std::vector<double> &phi , std::vector<double> &vel, int64_t &count, std::vector<int64_t> &start_idx, std::vector<int64_t> &end_idx, std::vector<int64_t> &pixnum_start, std::vector<int64_t> &pixnum_end, int rank_diff,  unordered_map<int64_t, int64_t>* ring_to_idx);

void get_pix_list_rank(int octant, int rank, int numranks, int64_t npix_lores, std::vector<int> pix_list_oct, std::vector<int> pix_list_noct,  std::vector<int> &lores_pix);

int compute_ranks_count( PLParticles* P, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int numranks, vector<int> &send_count, int rank_diff);

void clear_pixel(int64_t start_idx, int rank_diff, vector<double> &rho , vector<double> &phi, vector<double> &vel);

void compute_ranks_index( PLParticles* P, T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int numranks, vector<int> send_off, vector<int64_t> &id_particles, int rank_diff);

int assign_dm_ngp(vector<double> &rho, vector<double> &phi, vector<double> &vel, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx);

int assign_dm_cic(vector<double> &rho, vector<double> &phi, vector<double> &vel, PLParticles* P, T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx, float samplerate);


