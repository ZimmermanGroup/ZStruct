#ifndef RTYPE_H
#define RTYPE_H

#include "stringtools.h"
#include "utils.h"
#include "pTable.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>

using namespace std;

class RTYPE {

  private:

   int natoms;
   string* anames;
   int* anumbers;
   int* coordn;
   int nadd;
   int nbrks;
   int* add;
   int* brks;
   double Pr;

   int inited;
   int alloced;

   void sort_rxn();  
   void sort_rxn1(int nadd1, int* add1, int nbrks1, int* brks1, int* anumbers1, int* coordn1);


  public:

   int id;

   void init();
   void freemem();
   void set_rxn(int natoms1, string* anames1, int* anumbers1, int* coordn1, int nadd1, int* add1, int nbrks1, int* brks1);
   void set_pr(double Pr1);
   double get_pr();
   int get_nadd();
   int get_nbrks();

   int match(int natoms1, int* anumbers1, int* coordn1, int nadd1, int* add1, int nbrks1, int* brks1);
   void print();
 
};

#endif
