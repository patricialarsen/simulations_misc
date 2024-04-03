#include "healpix_base.h"
#include <unordered_map>
#include <string>
#include <vector>
using namespace std;

void get_pix_list_rank(int octant, int rank, int numranks, int64_t npix_lores, std::vector<int> pix_list_oct, std::vector<int> pix_list_noct,  std::vector<int> &lores_pix);

int compute_ranks_count_halos( PLHalos* P, T_Healpix_Base<int> map_lores,  int numranks, vector<int> &send_count, int rank_diff);


void compute_ranks_index_halo( PLHalos* P, T_Healpix_Base<int> map_lores, int numranks, vector<int> send_off, vector<int64_t> &id_particles, int rank_diff);

void clear_pixel(int64_t start_idx, int rank_diff, vector<float> &rho , vector<float> &phi, vector<float> &vel);
void clear_pixel_hydro(int64_t start_idx, int rank_diff, vector<float> &rho , vector<float> &phi, vector<float> &vel, vector<float> &ksz, vector<double> &tsz, vector<double> &xray1, vector<double> &xray2);


int assign_sz_cic_halocut(vector<float> &rho, vector<float> &phi, vector<float> &vel,vector<float> &ksz, vector<double> &tsz, PLParticles* P, PLHalos* H,  T_Healpix_Base<int64_t> map_hires, vector<int64_t> pixnum_start, vector<int64_t> pixnum_end, vector<int64_t> start_idx,  unordered_map<int64_t, int64_t> ring_to_idx, float hval, bool borgcube, bool adiabatic, float samplerate, string cloudypath, float masscut, float radiuscut);


