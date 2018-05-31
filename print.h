void print_xyz_gen(int natoms, string* anames, double* coords);
void print_triple_xyz(int natoms, string* anames, int* anumbers, double* c0, double* c1, double* c2, double* e);
void print_single_xyz_save(string xyzfile_string, int natoms, string* anames, int* anumbers, double* coords0);
void print_double_xyz_save(string file, int natoms, string* anames, double* c0, double* c1, double* e);
void print_triple_xyz_save(string file, int natoms, string* anames, double* c0, double* c1, double* c2, double* e);

void save_gephi(int npts1, int* ids, double* distances, double* values, double alpha, double savethresh);
void save_gephi_2(int npts1, int* ids, double* distances, double* values, double* values_print, double alpha, double savethresh);
void save_gephi_3(int npts1, int* ids, double* distances, double* values, double* values_print, double alpha, double savethresh);
