#ifndef RXNFTR_H
#define RXNFTR_H

#include "atom.h"
#include "icoord.h"
#include "stringtools.h"
#include "utils.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>

class RXNFTR {

  private:

   int nadd;
   int nbrk;
   int nchg;
   int nswp;

   ATOM* atoms;
   int new_atom;
   int max_atoms;

   int maxat;
   ATOM at1;
   ATOM at2;
   ATOM* at1d;
   ATOM* at2d;

   double atom_match(int e1, int c1, double q1, int add1, int brk1, int swp1, int nat1, ATOM* at1, ATOM at2);
   double atom_match(int e1, int c1, double q1, int add1, int brk1, int swp1, int nat1, ATOM* at1, ATOM at2, int* matched);
   void get_attached(ICoord ic1, int a1, ATOM& atd1);
   int get_attached_all(ICoord ic1, int a1, ATOM* atd);

  public:

   int p;
   int natoms;
   int quiet;

   void init(int nelem1);
   void reset(int nelem1);
   void update(int nelem1);
   void print_atoms();
   int count_parameters();

   double get_value_0(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1);
   double get_value_1(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1);

   void set_element(int wa, int e1);
   void set_coordn(int wa, int c1);
   void set_q(int wa, double q1);
   void set_qsig(int wa, double qsig1);
   void set_abs(int wa, int add1, int brk1, int swp1);
   int get_element(int wa);
   int get_coordn(int wa);
   double get_q(int wa);
   double get_qsig(int wa);
   void get_abs(int wa, int& add1, int& brk1, int& swp1);

   void set_attached(int wa, ATOM at1);
   void set_attached(int wa, int e1);
   void set_attached(int wa, int e1, int c1);
   void set_attached(int wa, int e1, int c1, double q1);
   void set_attaching(int wa, ATOM at2);
   void set_attaching(int wa, int e1);
   void set_attaching(int wa, int e1, int c1);
   void set_attaching(int wa, int e1, int c1, double q1);
   void delete_attached(int wa);
   void delete_attaching(int wa);
   int get_nattached(int wa);
   int get_nattaching(int wa);
   void get_attached(int wa, ATOM& at1);
   void get_attaching(int wa, ATOM& at2);

   void freemem();

};

#endif
