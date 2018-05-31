#ifndef KNNR_H
#define KNNR_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <math.h>
#include <stdio.h>
using namespace std;

#include "stringtools.h"
#include "utils.h"
#include "print.h"

class KNNR {
  
  private:
  
   double* distances;
   double* values;
   double* values_print;

   int fsize;
   int ptskip; //for getting nn
   int nlw;

   void get_distances(int npts, double* X);
   void filter_positive(int k, double* knnd, double* knnw, int* knn);

   int find_knn(double* X1, int k, int* knn, double* knnd);
   int find_knn(int pt, int k, int* knn, double* knnd);
   int check_unique(int pt);

   void determine_active();
   void reassign_mem(int npts1);

  public:

   double LOW_WEIGHT;
   double pthresh;
   double sumw;
   double sumd;
   int npts;
   int knn_N;
   int quiet;
   double* active; //active features
   int* udata; //use these data points

   double* X;

   int* ids;
   double* errlist;
   int* knnlist;

   double test_points(int k);
   double predict_point(int pt, int k);
   double predict_point(double* X1, int k);
   void load_values(int npts1, int fsize1, double* X1, double* y1);
   void load_values_print(int npts1, double* y1);

   double get_distance(double* X1, double* X2);

   void init();
   void reset_active();
   void freemem();


};

#endif
