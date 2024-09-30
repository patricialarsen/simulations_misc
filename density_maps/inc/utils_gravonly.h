#include "healpix_base.h"
#include <string>
#include <vector>
using namespace std;


void read_and_redistribute(string file_name, int numranks, PLParticles* P,  T_Healpix_Base<int> map_lores, T_Healpix_Base<int64_t> map_hires, int rank_diff, bool output_downsampled = false, float downsampling_rate = 1.0, string file_name_output = "");

void write_files(string outfile, string stepnumber,vector<int64_t> start_idx, vector<int64_t> end_idx, vector<int64_t> pix_nums_start, vector<double> rho, vector<double> phi, vector<double> vel, int64_t npix_hires);



