#define POSVEL_T float
#define ID_T int64_t
#include <vector>

using namespace std;

int compute_rank_fullsky(int pix_val, int numranks);

int compute_ranks_halo( Halos_test* H_1, T_Healpix_Base<int> map_lores, int numranks, vector<int> &send_count, vector<halo_properties>* halo_send);

int redistribute_halos(int numranks, vector<halo_properties> halo_send, Halos_test* H_1, vector<int> send_count);

int read_and_redistribute(string file_name, T_Healpix_Base<int> map_lores, int numranks, Halos_test* H_1);

int write_files(Halos_test &H_1, int l, vector<int> pixel_nums_rank, vector<int64_t> pixel_counts, string zrange_str, string step_s);

int write_files_slab(Halos_test &H_1,int rank,  vector<int> pixel_nums_rank, vector<int64_t> pixel_counts, string zrange_str, string step_s);

void read_halos(Halos_test *H0, string file_name);
