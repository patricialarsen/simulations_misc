#include "healpix_base.h"
#include <string>
#include <vector>
using namespace std;

void read_and_redistribute(string file_name, int numranks, Particles* P,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff);
void write_files(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<float> rho, vector<float> phi, vector<float> vel);

void write_files_hydro(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<float> rho, vector<float> phi, vector<float> vel, vector<float> ksz, vector<float> tsz);