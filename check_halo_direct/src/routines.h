bool comp_by_halo_dest(const halo_properties_test &a, const halo_properties_test &b);
bool comp_by_halo_tag(const halo_properties_test &a, const halo_properties_test &b) ;
bool comp_by_fof_mass(const halo_properties_test &a, const halo_properties_test &b) ;
bool comp_by_fof_x(const halo_properties_test &a, const halo_properties_test &b) ;
bool comp_by_fof_y(const halo_properties_test &a, const halo_properties_test &b) ;

bool comp_by_fof_z(const halo_properties_test &a, const halo_properties_test &b) ;
bool check_comp_halo(const halo_properties_test &a, const halo_properties_test &b);
int vec_to_rank(float x, float y, float z, float box_size);

bool comp_by_sod_dest(const sod_binproperties_test &a, const sod_binproperties_test &b) ;
bool comp_by_sod_id(const sod_binproperties_test &a, const sod_binproperties_test &b) ;
bool comp_by_sod_bin(const sod_binproperties_test &a, const sod_binproperties_test &b) ;
bool comp_by_part_dest(const particles_test &a, const particles_test &b) ;
bool comp_by_part_tag(const particles_test &a, const particles_test &b) ;
bool comp_by_part_id(const particles_test &a, const particles_test &b) ;

	

int compute_mean_float_dist(vector<float> *val1, vector<float> *val2, int num_halos , int num_halos2, string var_name, string var_name2, float lim);
int compute_mean_float(vector<float> *val1, vector<float> *val2, int num_halos , string var_name, string var_name2, float lim);
int  compute_mean_std_dist_halo(Halos_test H_1 , Halos_test H_2, float lim );
int compute_mean_std_dist_halos2(Halos_test H_1 , Halos_test H_2, float lim );
int  compute_mean_std_dist_sod(SODBins_test H_1 , SODBins_test H_2, float lim );


void read_halos(Halos_test &H0, string file_name, int file_opt) ;
void read_sodfiles(SODBins_test &H0, string file_name);
void read_particles(Particles_test &H0, string file_name);


int perform_halo_check(string fof_file, string fof_file2, float lim, float min_mass, float max_mass, map<int64_t,int> *tag_map);
int match_pos (string fof_file, string fof_file2, float lim, float box_size, float min_mass, float max_mass);
int compare_dist(string fof_file,string fof_file2, float lim);
int sodbin_check(string fof_file, string fof_file2, float lim, map<int64_t,int> *tag_map);
int part_check(string fof_file, string fof_file2);

