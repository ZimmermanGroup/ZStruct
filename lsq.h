#ifndef LSQ_H
#define LSQ_H

#include "utils.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>

class LSQ {

  private:

   int p; //dimensionality
   int N; //number of data points


  public:

   int quiet;

   double* params;

   void go_test();

   void init();
   double do_least_squares(double* x, int p, double* y, int N);
   double estimate_point(double* x);

};

#endif
