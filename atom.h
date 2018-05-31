#ifndef ATOM_H
#define ATOM_H

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

class ATOM {

  private:

    int inited;

    int add;
    int brk;
    int swp;

    int elem;
    int coordn;

    double q; 
    double qsig;

    int nattached;
    int nattaching;
    ATOM* attached;
    ATOM* attaching;

    int element_match(int e1);
    int coordn_match(int c1);
    int add_match(int add1);
    int brk_match(int brk1);
    int swp_match(int swp1);
    int attached_match(ATOM at1);
    double attaching_match(ATOM at2);
    double q_eval(double q1);

  public:

   int p;
   int active;

   void init();
   void reset();
   void set_element(int e1);
   void set_coordn(int c1);
   void set_q(double q1);
   void set_qsig(double qsig1);
   void set_abs(int add1, int brk1, int swp1);
   void set_attached(ATOM at1);
   void set_attaching(ATOM at2);
   void set_attached(int e1);
   void set_attaching(int e1);
   void set_attached(int e1, int c1);
   void set_attaching(int e1, int c1);
   void set_attached(int e1, int c1, double q1);
   void set_attaching(int e1, int c1, double q1);
   void delete_attached();
   void delete_attaching();
   int get_element();
   int get_coordn();
   double get_q();
   double get_qsig();
   void get_abs(int& add1, int& brk1, int& swp1);
   int get_nattached();
   int get_nattaching();
   void get_attached(ATOM& at1);
   void get_attaching(ATOM& at2);


   int compare(int e1, int c1, int nadd1, int nbrk1, int nswp1);
   double compare(int e1, int c1, double q1, int nadd1, int nbrk1, int nswp1);
   double compare(int e1, int c1, double q1, int nadd1, int nbrk1, int nswp1, int nat1, ATOM* at1, ATOM at2);
   int count_parameters();
   void print();

   void freemem();

};

#endif

