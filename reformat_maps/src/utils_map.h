#define POSVEL_T float
#define ID_T int64_t
#include <vector>

using namespace std;

//int compute_ranks( MapData* M, int numranks, vector<int> &send_count, vector<map_properties>* maps_send);
int compute_ranks( MapData* M, int numranks, vector<int> &send_count, vector<map_properties> &maps_send);
//int redistribute_maps(int numranks, vector<map_properties> maps_send, vector<map_properties>* recv_maps, MapData* M, vector<int> send_count);
int redistribute_maps(int numranks, vector<map_properties> maps_send, vector<map_properties>& recv_maps, MapData* M, vector<int> send_count);

//int redistribute_maps(int numranks, vector<map_properties> maps_send, MapData* M, vector<int> send_count);

int read_and_redistribute(string file_name, int numranks, MapData* M, vector<map_properties>& maps_send, vector<map_properties>& maps_recv);

int write_files_slab(MapData &M, int rank, string step_s, vector<float>& out_float);

void read_map_file(MapData* M, string file_name);

//void read_map_file(MapData* M, string file_name);

