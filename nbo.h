#ifndef NBO_H
#define NBO_H

#include "stringtools.h"
#include "utils.h"

class NBO {
  
  private:
  
   int natoms;
   int* anumbers;
   string* anames;
   double* xyz;

   int isQC; //QC or MOPAC
   int mheadersize;
   string* mheader;

   int nalpha; 
   int nbeta;
   int nelec;
   int nao; //# atomic orbitals
   int nmo; //# molecular orbitals
   int read_mo(string filename);
   int read_mopac_mo(string filename);
   void clean_long_bonds();
   double distance(int i, int j);
   void get_three_center(int& cmo, vector<string> tok_line);
   void alloc_nbo();

  public:

   void alloc(int natoms);
   void init(int natoms0, int* anumbers0, string* anames0, double* xyz0);
   void reset(int natoms0, int* anumbers0, string* anames0, double* xyz0);
   void freemem();

   int read_nbo_file(string filename);
   void print_nbo();

   double* q;
   double* spratio;
   double* MO;
   int* wAO;
   int* tAO;

   int nb;
   int* blist; //broken list

   int hasNBO;
   int bmo; //occupied bonding orbitals
   int vmo; //unoccupied bonding orbitals
   double* mo_occ;
   double* bmo_occ;
   double* bmo_polar;
   int* bmo_atoms;
   int* bmo_num;
   double* vmo_occ;
   int* vmo_atoms;
   int* vmo_num;
   
   void print_molden_orbs(int norbs, int* olist, string filename);
   int compare_nbo(NBO nbo1, string filename, int quiet);
   double get_pol(int a1, int a2, double& ev1);

};

#endif
