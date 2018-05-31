#ifndef DFT_H
#define DFT_H

#include "stringtools.h"
#include "pTable.h"

class DFT {
  
  private:
  
   int natoms;
   int* anumbers;
   string* anames;

   void xyz_read(string filename);

  public:

   double sp();
   double sp(string filename);
   void sp_dnr(string filename);
   void opt_dnr(string filename);
   double get_energy(string filename);
   double get_opt_energy(string filename);
   void get_structure(string filename, double* xyz);
   void get_charges(string filename, double* q);

   double ts();
   double ts(string filename);
   void ts_dnr(string filename);
   double get_energy_ts(string filename);

//   double opt();
//   double opt(string filename);
   void alloc(int natoms);
   void init(int natoms, int* anumbers, string* anames, double* xyz);
   void reset(int natoms, int* anumbers, string* anames, double* xyz);
   void freemem();

   double energy0;
   double energy;
   double energyts;
   int converged;

   double* xyz0;
   double* xyz;

};

#endif
