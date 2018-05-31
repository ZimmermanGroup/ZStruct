#ifndef ZSTRUCT_H
#define ZSTRUCT_H

#define MAX_ADD 2
#define MAX_BRK 2
#define MAX_CHG 3
#define MIN_CHG 1
#define RINGSHUTTLE 4
#define USE_H_TRANSFER 0
#define DO_PAIR_ISOS 1
#define TM_PAIRS 0

#define DO_NOT_WRITE 0
#define DO_NOT_WRITE_REACTS 1
#define GSM_FRAC_DONE 0.99
#define USE_MOPAC 0

#define USE_NBO 1
//CPMZ check
#define ISOMER_DB 0
#define USE_DT 0
#define USE_GA 0
#define GA_POP 16
#define GA_STEPS 1
#define TEST_GA 0
#define ALL_PASS 1

//standard library includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <string>
#include <cctype>
#include <ctime>
#include <sys/stat.h>

#include "stringtools.h"
#include "icoord.h"
#include "mopac.h"
#include "rxndb.h"
#include "align.h"
#include "dft.h"
#include "print.h"
#include "utils.h"
#include "nbo.h"
#include <omp.h>


using namespace std;


struct irxn1
{
  int id; //internal reference #

  int natoms;
  int nadd;
  int nbrks;
  int* adds;
  int* brks;

  double* brkd;

//  string* addelm;
//  string* brkelm;

  int* anums;
  int* coordn;
  string* anames;
  double* xyz;

  double estart;
  double ets;
  double erxn;

  double* q;
};


class ZStruct
{

 private:

  int nreact;
  int nreacta;
  int* reactnum;
  int* reactrun;
  int npair;
  int* pairs; 
  int* pairsnum;
  int* rsaved;
  int wpaira;
  int* wpair; //over niso
  int* wshpair;
  int* upair;
  NBO* nbo_reacts;
  int ga_avail;
  RXNDB rxns1;
  RXNDB rxns2;
  double pthresh;
  int atomlimit;

  int niso;
  int nisoa; //# alloc'd

  int nelem;
  int* climit_l;
  int* climit_h;
  int activea;
  int** active;
  int nhh;

  int* pair_react;
  int* natomsr;
  double** qr; //over nreact
  double** qp; //over niso

  int* natoms;
  string** anames;
  int** anumbers;
  double** xyzr;
  double** xyzp;
  double** xyzt;
  double** qrc;

  int elista;
  double* elistref; //reference energies for Ea
  double* elistr; //individual reactant energies
  double* elistrp; //separated pair energies
  double* elistc; //complexation energies
  double* elistcmin; //minimum of complexation energies for pairs
  double* elistcminsh; //including shuttles
  double* elistts; //ts energies
  double* elistp; //product energies
  double* Ea; //activation energies
  double* Erxn; //energy of reaction

  Align align1;
  int nfragb;
  int* fragb;

  int nshuttle;
  int* wshuttle;
  int** shuttles;
  ICoord* icshuttle;

  void init();
  void dealloc();

  int read_reactants(ICoord* icr);
  int pair_allowed(int i, int j, ICoord* icr);
  int get_frozen(int nreact);
  int get_limits();
  int hh_bond(int a1, int a2, int* anumbers, int* coordn);
  int two_frags(int a1, int a2, ICoord ic1);
  int two_frags(int a1, int a2, int a3, int a4, ICoord ic1);
  double get_tse(int nnodes, int& nnmax, double* energies);
  void create_frag_bonds(ICoord ic1);

  int read_shuttles();
  void save_shuttles(int type);
  void read_shuttles_data();
  int is_shuttle(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1, int* s1);
  int get_shuttle(ICoord ic1, int nr, int wsh, int nadd1, int* add1, int nbrks1, int* brks1, int* s1);

  void fill_elistrp();
  void update_elist_mem(int start, int end);
  void update_react_mem();
  void update_pairs_mem(int npairs0, int npairp);
  void update_wpair_mem(int nisop, int nisof);
  void update_active_mem(int nreact1);
  void get_elistcmin(int wp, int niso1);
  void get_elistref();
  int formula_match(int ws1, int ws2);
  void get_pair_indices(int wp, int niso1, int& pstart, int& pend);
  void assign_sequential_pairs(int wp);

  void write_duplicates(int nsave);
  int read_duplicates();
  int read_gsm_data();
  void save_gsm(int type);
  void save_frozen(int natoms1, int* active1, int num);
  double read_temperature();
  int read_string(int wfile, double* energies, int natoms0, double** xyz);
  void read_gsm_files(int first, int last, int* gsmdone);
  void write_initial_xyz(int wfile, int natoms, string* anames, double* xyz, int q1);
  void write_pair_xyz(int wfile, int natoms, string* anames, double* xyz);
  void write_ISOMER(int wfile, int nadd, int* add, int nbrk, int* brk);
  void write_ISOMER(int wfile, int nadd, int* add, int nbrk, int* brk, int nangle, int** angles, double* anglev, int ntor, int** torsion, double* torv);
  int read_ISOMER(int wfile, int& nadd, int* add, int& nbrks, int* brks);

  void dft_para(int nreact, ICoord* icr);
  void dft_para_qchem(int nreact, ICoord* icr);
  void dft_para_ase_gaussian(int nreact, ICoord* icr);
  void write_charge(int q1, string nstr);
  void dft_para_q(int niso1);
  void get_nbo_reactants(ICoord* icr);
  void get_nbo_reactions(ICoord* icr);
  void get_one_nbo(int a1, int type, NBO* ntest1, Mopac& mopt1);
  void gsm_para(int first, int last);
  void gsm_analyze(int wp, int first, int last, string filename);
  void get_qrc(int ws, double* qrc1);
  int analyze_path(int ws, int wp, int nadd, int* add, int nbrks, int* brks, int nnodes, double* energies, double** xyz, ICoord ic1);
  void store_xyz(int ws, int nnodes, int natoms0, string* anames0, int* anumbers0, double** xyz0);

  void form_rxn_ga();
  void write_ga();
  int read_ga();

  int create_react_isos(ICoord* icr);
  int create_pairs_isos(ICoord* icr);
  int create_isos(ICoord ic1, int nr, int niso1, int* active1, int wsh);
  int get_breaks(int nbrks, int* brks, ICoord ic1, int* active1);
  int get_combo(int nadd, int nbf, int nbrks, int* brks, ICoord ic1, int nr, int* active1, int wsh);
  int get_combo_h1(int nadd, int nbf, int nbrks, int* brks, ICoord ic1, ICoord ic1tm, int nr, int* active1, int wsh);
  int get_h_transfer(ICoord ic1, int nadd, int* add, int nbrks, int* brksp, int natoms, int* coordn, int* moving, int wsh);

  int tm_angle_drive(int nr, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* active1);
  int determine_geomtype(int wtm, int stm, int coordntm, ICoord ic1);
  int run_tm(int coordntm, int geomtype, int wtm, int stm, int naddtm, int* addtm, int nbrktm, int* brktm, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1);
  int run_2c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv);
  int run_2c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv);
  int run_3c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv);
  int run_3c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv);
  int run_4c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1);
  int run_4c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1);
  int run_5c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1);
  int run_5c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1);
  int run_6c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1);

 //organometallic build
  int run_4c_om(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv);

  int find_2c(int wtm, int stm, ICoord ic1);
  int find_3c(int wtm, int stm, ICoord ic1);
  int find_4c(int wtm, int stm, ICoord ic1);
  int find_5c(int wtm, int stm, ICoord ic1);
  int find_6c(int wtm, int stm, ICoord ic1);
  int find_opp(int a1, int wtm, int stm, ICoord ic1);
  void find_linear_pair(int& a1o, int& a2o, int& a1f, int& a2f, int a1, int a2, int a3, int a4, int wtm, ICoord ic1);
  void find_2_at_90(int& a1s1, int& a1s2, int a1, int a1o, int wtm, int stm, ICoord ic1);
  int find_on_top_5c(int wtm, int stm, ICoord ic1);
  void find_tb_5c(int* a, int wtm, int stm, ICoord ic1);
  void find_oct_planes(int wtm, ICoord ic1, int& p11, int& p12, int& p13, int& p14, int& p21, int& p22, int& p23, int& p24, int& p31, int& p32, int& p33, int& p34);
  void find_5c_angles(int nattached, int* attached, int wtm, ICoord ic1, double* angles);
  void find_6c_angles(int nattached, int* attached, int wtm, ICoord ic1, double* angles);
  void find_linear_pairs(int& a1o, int& a2o, int& a3o, int& a4o, int a1, int a2, int a3, int a4, int wtm, ICoord ic1);
  void get_perp_vec(int a1, int a2, double* coords, double* w);
  void get_oop_vec(int a1, int a2, int a3, double* coords, double* w);
  void get_vec(int a1, int a2, double* coords, double* w);
  void get_Y_vec(int a1, int a2, int a3, double* coords, double* w);
  void get_vec_45p_align(int a1, int a2, int a3, int a4, double* coords, double* w);
  void get_vec_45m_align(int a1, int a2, int a3, int a4, double* coords, double* w);
  void get_Y_vec_45p(int a1, int a2, int a3, double* coords, double* w);
  void get_Y_vec_45m(int a1, int a2, int a3, double* coords, double* w);
  void get_vec_45p(int a1, int a2, int a3, double* coords, double* w);
  void get_vec_45m(int a1, int a2, int a3, double* coords, double* w);
  void align_vec(int a1, int a2, double* coords, double* w);
  void sort_angles(int nangles, double* angles);
  int eta2_double_break(int wtm, int stm, int ncoord, int nbrks, int* brks, ICoord ic1);
  int find_attached(int wat, int stm, ICoord ic1, int* attached);
  int get_wtm_addbrktm(int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int& naddtm, int* addtm, int& nbrktm, int* brktm);


  int ml_eval(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1);
  int ga_eval(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1);

  int diff_structure(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  int diff_structureiq(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  int diff_structureigq(int natoms1, int natoms2, string* anames1, string* anames2, int* anumbers1, int* anumbers2, double* xyz1, double* xyz2);
  int diff_structureq(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  int diff_structurec(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  int diff_structurecq(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  void swap_atoms(double* xyz, int a1, int a2);

  void do_one_step(ICoord* icr);
  void create_new_reactants(int* rxncont, ICoord* icr);
  void calc_Ea_Erxn();
  void form_reaction_network(int* rxncont, ICoord* icr);
  void mark_duplicates(ICoord* icr);

  int check_against_reactants(int natoms1, string* anames1, int* anumbers1, double* xyz1, ICoord* icr);
  void save_new_reactant(ICoord* icr, int natoms1, string* anames1, int* anumbers1, double* xyz1, int* active1);

//CPMZ full DB construction
  void store_isomers(int id1, int wp1, int nadd1, int nbrks1, int* add1, int* brks1, double* E1, int natoms1, string* anames1, int* anums1, double* xyz1, double* qrc1, int* coordn1);
  void save_allrxns();
  vector<irxn1> allrxns;

 public:

  double temperature;
  void go_zstruct(int nsteps);

};


#endif

