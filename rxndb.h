#ifndef RXNDB_H
#define RXNDB_H

#include "rxnftr.h"
#include "icoord.h"
#include "stringtools.h"
#include "utils.h"
#include "lsq.h"
#include "knnr.h"
#include "rtype.h"
#include "nbo.h"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>

class RXNDB {

  private:

   int level; //analysis technique

   int N;
   int* p;
   int plus_b;
   int Nalloc; //# allocated
   int ng;
   int maxg;
   int nbest_fv;
   int* best_fv;
   int nrxnftrs;
   RXNFTR** rxnftrs;
   double* ga_error;
   int ga_inited;
   NBO* nbosr;
   NBO* nbosp;

  //data for each rxn
   int* natoms;
   string** anames;
   int** anumbers;
   double** xyzr;
   double** xyzp;
   int** coordnr;
   int* nadd;
   int* nbrks;
   int** add;
   int** brks;
   double* Ea;
   double* Erxn;
   double** qrc; //reactant complex charges from individual reactants
   double* Pr;
   double* PrT2;
   int lastAdd;

   int nelem; //max elements
   int nelemf; //# found
   int* emap;
   int* remap;
   int* efindex;
   int* climit_l;
   int* climit_h;

   int naddbrktypes;
   double* allow_add_brk_elem;
   double* addbrk_layer;

   int create_features_from_db(int wg);
   int create_features_climit(int wg);
   int create_features_climit_ab(int wg);
   int create_features_climit_abs(int wg);
   int create_features_climit_abs_att(int wg);
   int create_features_climit_abs_att_2(int wg);
   int count_elements();
   int count_parameters(int wg);
   void create_sigmoid();
   void get_rowX_nbo(int wg, int wr, double* X, int* M, int sort);
   void get_rowX(int wg, int wr, double* X);
   void get_rowX(int wg, double* X, ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1);
   double get_rowX_norm(int wg, double* X, ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1);

   void selection(int nbestg, int* bestg);
   void spawn_and_splice(int nbestg, int* bestg);
   int mutate(int wg);
   void mutate_one(int wg, int r1);
   void copy_atom(RXNFTR& ftr1, RXNFTR ftr0);
   int replace_dead_features(int wg);

   int get_abtype(int nadd1, int nbrks1);
   int get_abtype(int i, int nadd1, int nbrks1);

   void alloc_db(int natoms);
   void write_xyz(int wfile, int wtype, int natoms1, string* anames1, double* xyz1, double Ea1);
   int check_duplicate(int ws0, int natoms0, int* anumbers0, double* energies,
           int nadd0, int* add0, int nbrks0, int* brks0);
   void reorder_nbo_db();
   void fix_unbroken_bonds();
   void fix_unadded_bonds();

   LSQ* lsqs;
   KNNR* knnrs; 
   int* nrtypes;
   RTYPE** trxns;


  public:

   int quiet;
   double temperature;
   int* ids;
   double* pthresh;
   double* similar; //when using add_ts_xyz_test
   int* unique; 
 
   void init(int level0);
   void set_limits(int* climit_l, int* climit_h);
   void freemem();

   void init_ga(int ng);
   void init_ga_2();
   void init_ga_features(int ng);
   void copy_features(int wg1, int wg0);
   void copy_features(RXNFTR* &rxnftrs1, RXNFTR* rxnftrs0, int size1, int size0);
   void ga_run(int ncycles);

   void write_ga(string filename);
   int read_ga(string filename);

   int add_ts(int ws, int natoms0, string* anames0, int* anumbers0, double* xyzr0, double* xyzp0, double* xyzts0, double* energies,
           int nadd0, int* add0, int nbrk0, int* brks0, double* qrc0);
   int add_ts_xyz(int ws, int natoms0, string* anames0, int* anumbers0, double* xyzr0, double* xyzp0, double* xyzts0, double* energies, double* qrc0);
   int add_ts_xyz_test(int ws, int natoms0, string* anames0, int* anumbers0, double* xyzr0, double* xyzp0, double* xyzts0, double* energies, double* qrc0,
           int nadd0, int* add0, int nbrk0, int* brks0);
   int find_ts(int ws);
   double create_regression(int wg);
   double create_regression_nbo(int pr);
   void find_pthresh(int wg, double pfalseneg);
  
   void add_nbo_data(int nnbo, NBO* nbo1r, NBO* nbo1p, int* ids1);
   void save_xy_data(string filename);
   int read_xy_data(double* X1, double* Ea1, double* Erxn1, int* ids1, string filename);
   void set_unique(int N1, int* unique1, int* ids1);

   void element_analysis(int wg);
   void make_decision_tree();
   void make_combo_tree(int naddmax, int nbrkmax);

   int tree_screen_0(int nadd0, int nbrks0);
   double tree_screen_1(int natoms0, string* anames, int* anumbers, double* xyzr, int nadd0, int* add0, int nbrks0, int* brks0, double* qrc0);
   double tree_screen_2(int natoms0, string* anames, int* anumbers, int* coordn, double* xyzr, int nadd0, int* add0, int nbrks0, int* brks0, double* qrc0);
   double eval_prob(int natoms0, string* anames, int* anumbers, double* xyzr, int nadd0, int* add0, int nbrks0, int* brks0, double* qrc0);
   void print_reactions();
   void print_reactions_2();

   void save_structures();

};

#endif
