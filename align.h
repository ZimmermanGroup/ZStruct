#ifndef ALIGN_H
#define ALIGN_H

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

void print_xyz_gen(int natoms, string* anames, double* coords);

class Align {

  private:

//   int nadd;
//   int* add;

   string* anames1;
   string* anames2;
   int* anumbers1;
   int* anumbers2;
   double* xyz1s;
   double* xyz2s;
   double* xyz1;
   double* xyz2;

   int natoms3;
   string* anames3;
   int* anumbers3;
   double* xyz3;

   int check_frag(int a1, int a2);
   int check_frag_3(int a1, int a2);
   int get_bonds(int a1, ICoord ic1, int* bonded);
   void align_to_x(int natoms, int t1, int t2, double* xyz, string* anames, int sign, double offset);
   void rotate_around_x(int natoms1, int natoms2, double torv, double* xyz);
   void linear_right(double* v1, int a1, int* bonded, double* xyz);
   void planar_cross(double* v1, int a1, int* bonded, double* xyz);
   void align_v1(int nvf, double* v1);
   void point_out(double* v1, int natoms, double* xyz);
   void vdw_vector_opt(int natoms1, int natoms2, double* v1, ICoord icp);


  public:

   int inited;
   double* xyza;
   double* xyza3;
   int natoms1;
   int natoms2;
   int avec1; //final rotation alignment, not implemented
   int avec2;

   void init(int natoms1i, string* anames1i, int* anumbers1i, double* xyz1i, int natoms2i, string* anames2i, int* anumbers2i, double* xyz2i);
   void add_third(int natoms3i, string* anames3i, int* anumbers3i, double* xyz3i);
   void align_zero();
   void add_align(int nadd1, int* add1);
   int add_align_v(int nadd1, int* add1, int wtm, double* aprv);
   void shuttle_align(int nadd1, int* add1);

   void print_xyz();
   void print_xyz_3();

   void freemem();

};

#endif

