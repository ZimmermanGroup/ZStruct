#ifndef ICOORD_H
#define ICOORD_H

#include "stringtools.h"
#include "pTable.h"

class ICoord {

  private:

    int data;
//    int* atypes;                  //array of MM atom types  // assign values directly?
    double* amasses;              //array of atomic masses
//    double* charge;               //array of MM atomic charges
    double* dxm1;

    double* bondd;

    double* anglev;

    double* torv;

    double* imptorv;
    int** imptor;

    int max_bonds;
    int max_angles;
    int max_torsions;
    int max_imptor;

    int max_nonbond;
    int n_nonbond;
    int** nonbond;
    double* nonbondd;

  void structure_read(string xyzfile);
  void alloc_mem();
  void make_bonds();
  void split_h();
  void coord_num();
  void make_angles();
  void make_torsions();

  void get_eta(int& tbond, int** btm);
  void make_imptor(); 
// ic creation function when bonds are already made
  void make_imptor_nobonds(); 

  int make_nonbond();

//
  void update_bonds();
  void update_angles();
  void update_torsion();
  void update_imptor();
  void update_nonbond();

  void create_xyz();

  // MM force field params
  double* ffR;
  double* ffeps;
  double ffbondd(int i, int j);
  double ffbonde(int i, int j);
  double ffangled(int i, int j);
  double ffanglee(int i, int j);
  double fftord(int i, int j, int k, int l); 
  double fftore(int i, int j, int k, int l); 
  double fftorm(int i, int j, int k, int l); //multiplicity 
  double ffimptore(int i, int j, int k, int l); 
  double ffimptord(int i, int j, int k, int l); 
  // function to make arrays?

  // Gradient terms
  double* grad;
  double pgradrms;
  void print_grad();
  void bond_grad_all();
  void bond_grad_1(int i, int j);
  double bond_stretch(int i, int j);
  void angle_grad_all();
  void angle_grad_1(int i, int j, int k);
  void torsion_grad_all();
  void torsion_grad_1(int i, int j, int k, int l);
  void imptor_grad_all();
  void imptor_grad_1(int i, int j, int k, int l);
  void vdw_grad_all();
  void vdw_grad_1(int i, int j);
  double vdw_energy_all();
  double vdw_energy_1(int i, int j);

  //Optimizer
  void update_xyz_sd();
  void update_xyz_cg();

//  ofstream xyzfile;
  void print_xyzf(ofstream xyzfile); // print xyz coords to file

  public:

  string comment;
  int ic_create_tm();
  int** bondstm;
  int nbondstm;
  int coordntm;
  int** bondstm2;
  int nbondstm2;
  int coordntm2;

  int** bonds;
  int nbonds;
  int** angles;
  int nangles;
  int** torsions;
  int ntor;
  int q1; //charge of this species
  int s1; //multiplicity
  double* q;

  double gradrms;
    
  int id; //for geoms[id] in zstruct
  int pid; // previous structure id
  double seenergy;
  double segsmenergy;
  double dftenergy;
  double dftlstenergy;
  double dftgsmenergy;
  double dfttsenergy;

  int natoms;
  double* coords;
  double* coordsts;
  double* coords0;
  string* anames;               //array of atomic symbols 
  int* anumbers;                //array of atomic indices 
  int* coordn;                  //coordination number
  int nimptor;

  int ic_create();
  int ic_create_nobonds();
  double mm_energy();
  int mm_grad();
  int mm_grad(ICoord shadow);
  int opt();
  int opt(string xyzfile, ICoord shadow);
  int opt(string xyzfile);
  void update_ic();
  void mm_init();

// help functions for iso
  int bond_exists(int b1, int b2);
  int bond_num(int b1, int b2);
  int hpair(int a1, int a2);
  int h2count();
  int isTM(int anum);

  int same_struct(double* xyz);

  int init(string xyzfile);
  int init(int natoms, string* anames, int* anumbers, double* xyz);
  int alloc(int size); 
  int reset(int natoms, string* anames, int* anumbers, double* xyz);
  int reset(double* xyz);
  void print_ic();
  void print_bonds();
  void print_xyz();
  void print_xyz_save(string filename);
  void print_xyz_save_no_natoms(string filename);
  void print_xyz_save(string xyzfile_string, double energy);

  int nfrags;
  int* frags;
  void make_frags();
  void bond_frags();

  double farBond;

  double getR(int i);
  double distance(int i, int j);
  double angle_val(int i, int j, int k);
  double angle_val_eta2(int i, int j, int k1, int k2);
  double angle_val_eta2_eta2(int i1, int i2, int j, int k1, int k2);
  double torsion_val(int i, int j, int k, int l); // for imptor and torsion

  void freemem();


};



#endif

