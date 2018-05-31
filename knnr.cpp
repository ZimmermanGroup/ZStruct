// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "knnr.h"

//was using 0.25 and 0.25/1.0
//note using 1.0 leads to many out of set hits
#define KNNR_ALPHA 0.15
#define CZERO 0.15
#define DZERO 0.05
#define TEST_MODE 1
#define PRINT_DISTS 1
#define POSITIVE_BIAS 0
#define TAKE_POSITIVE_MATCH 0


//Notes:
// 1. get_distances only operates over active features (from training set)


void KNNR::filter_positive(int k, double* knnd, double* knnw, int* knn)
{
  int wpt0 = knn[0];
  int wpt1 = knn[1];
  if (knnd[0]<CZERO)
  if (values[wpt0]>2.0*values[wpt1])
  {
    //printf(" fp knnd[0]: %4.3f \n",knnd[0]);
    knnw[0] = 1.0;
    for (int i=1;i<k;i++)
      knnw[i] = 0.;
  }

  return;
}

double KNNR::test_points(int k)
{
  if (k<1)
    return 99.;
  if (values==NULL)
  {
    printf(" in test_points, value not init'd \n");
    return 99.;
  }
  //printf(" in test_points, npts: %i \n",npts);

  double error = 0.;
#if TEST_MODE
  if (errlist!=NULL) delete [] errlist;
  errlist = new double[npts];
  if (knnlist!=NULL) delete [] knnlist;
  knnlist = new int[npts*k];
#if PRINT_DISTS
  get_distances(npts,X);
#endif
#endif

  nlw = 0;
  int nunique = 0;
  int npositive = 0;
  int ndistant = 0;
  double val;
  double err1 = 0.;
  for (int i=0;i<npts;i++)
  if (check_unique(i))
  {
    val = predict_point(i,k);
    if (val>pthresh) npositive++;
    if (val==1.2) ndistant++;
    err1 = val - values[i];

#if 1
    if (values[i]>pthresh)
    {
      if (val<pthresh)
        error += 2.0; //false neg penalty
      else
        error += 0.5* err1*err1;
    }
    else
    {
      error += err1*err1;
    }
#else
    error += err1*err1 * (1+1*values[i]);
#endif

#if TEST_MODE
    errlist[i] = err1;
#endif

    if (!quiet)
    {
      printf(" %4i actual pred: %6.3f %6.3f",i,values_print[i],val);
      printf("  %6.2f  %6.4f",sumd,sumw);
#if 0
      if (values[i]>0.2 && val < values[i]/2.5)
      {
        printf(" *(%6.5f)",sumd);
        if (ids!=NULL) printf(" %3i",ids[i]);
      }
#endif
      printf("\n");
    }

    nunique++;
  } //loop i over data points
  error = sqrt(error/npts);


  if (!quiet && active!=NULL)
  {
    printf(" active features:");
    for (int i=0;i<fsize;i++)
      printf(" %1.0f",active[i]);
    printf("\n");
  }
  int npositiver = 0;
  for (int i=0;i<npts;i++)
  if (check_unique(i) && values_print[i]>pthresh)
    npositiver++;
  double passr = 100.*npositiver/nunique;
  double pass = 100.*npositive/nunique;
  if (!quiet)
    printf(" %3i of %3i pass: %3.1f%%  true positive: %3.1f%% (%4i) ndistant: %4i \n",npositive,nunique,pass,passr,npositiver,ndistant);

  return error;
}




double KNNR::predict_point(double* X1, int k)
{
  double val = 0.;

  int* knn = new int[k];
  double* knnd = new double[k];
  for (int i=0;i<k;i++) knn[i] = -1;
  for (int i=0;i<k;i++) knnd[i] = 10000.;

  int nfound = find_knn(X1,k,knn,knnd);

#if TEST_MODE
  for (int i=0;i<k;i++)
    knnlist[i] = knn[i];
#endif

#if TAKE_POSITIVE_MATCH
  for (int i=0;i<k;i++)
  if (knnd[i]<CZERO)
  if (values[knn[i]]>pthresh)
  {
    val = values[knn[i]];

    delete [] knn;
    delete [] knnd;

    sumd = 0.;
    sumw = 5.; 
    return val;
  }
#endif

#if 0
  double totdist = 0.;
  for (int i=0;i<k;i++)
    totdist += knnd[i];
  for (int i=0;i<k;i++)
  if (knnd[i] < 0.00001)
    totdist = 0.;
#endif

  double* knnw = new double[k];

#if 1
  //exponential weighting
  double alpha = KNNR_ALPHA;
  for (int i=0;i<k;i++)
    knnw[i] = exp(-knnd[i]/alpha);
#else
 //inverse weighting  
  for (int i=0;i<k;i++)
  if (knnd[i] > 0.00001)
    knnw[i] = 1./knnd[i];
  else
    knnw[i] = 1000000.;
#endif

#if POSITIVE_BIAS
  filter_positive(k,knnd,knnw,knn);
#endif

  sumd = 0;
  for (int i=0;i<k;i++)
    sumd += knnd[i];
  sumw = 0.;
  for (int i=0;i<k;i++)
    sumw += knnw[i];
  for (int i=0;i<k;i++)
    knnw[i] = knnw[i] / sumw;

  if (sumw < LOW_WEIGHT)
  {
    //printf(" low weights \n");
    nlw++;

    delete [] knn;
    delete [] knnd;
    delete [] knnw;

    return 1.2;
  }

#if 0
  printf(" weights:");
  for (int i=0;i<k;i++)
    printf(" %4.3f",knnw[i]);
  printf("\n");
#endif

  for (int i=0;i<k;i++)
    val += knnw[i] * values[knn[i]];


  delete [] knn;
  delete [] knnd;
  delete [] knnw;

  return val;
}


double KNNR::predict_point(int pt, int k)
{
  double val = 0.;

  int* knn = new int[k];
  double* knnd = new double[k];
  for (int i=0;i<k;i++) knn[i] = -1;
  for (int i=0;i<k;i++) knnd[i] = 10000.;

  int nfound = find_knn(pt,k,knn,knnd); 

  for (int i=0;i<k;i++)
  if (knn[i] == pt)
  {
    printf(" ERROR: knn==pt: %i %i \n",knn[i],pt);
    exit(1);
  }

#if TEST_MODE
  for (int i=0;i<k;i++)
    knnlist[k*pt+i] = knn[i];
#endif

#if TAKE_POSITIVE_MATCH
  for (int i=0;i<k;i++)
  if (knnd[i]<CZERO)
  if (values[knn[i]]>pthresh)
  {
    val = values[knn[i]];

    delete [] knn;
    delete [] knnd;

    sumd = 0.;
    sumw = 5.;
    return val;
  }
#endif

#if 0
  double totdist = 0.;
  for (int i=0;i<k;i++)
    totdist += knnd[i];
  for (int i=0;i<k;i++)
  if (knnd[i] < 0.00001)
    totdist = 0.;
#endif

  double* knnw = new double[k];

#if 1
  //exponential weighting
  double alpha = KNNR_ALPHA;
  for (int i=0;i<k;i++)
    knnw[i] = exp(-knnd[i]/alpha);
#else
 //inverse weighting  
  for (int i=0;i<k;i++)
  if (knnd[i] > 0.00001)
    knnw[i] = 1./knnd[i];
  else
    knnw[i] = 1000000.;
#endif

#if POSITIVE_BIAS
  filter_positive(k,knnd,knnw,knn);
#endif

  sumd = 0;
  for (int i=0;i<k;i++)
    sumd += knnd[i];
  sumw = 0.;
  for (int i=0;i<k;i++)
    sumw += knnw[i];
  for (int i=0;i<k;i++)
    knnw[i] = knnw[i] / sumw;

#if 0
  if (pt==151)
  {
  printf("\n sumw: %4.3f \n",sumw);
  printf(" distances:");
  for (int i=0;i<k;i++)
    printf(" %4.3f",knnd[i]);
  printf("\n");
  printf(" weights:");
  for (int i=0;i<k;i++)
    printf(" %4.3f",knnw[i]);
  printf("\n");
  printf(" values:");
  for (int i=0;i<k;i++)
    printf(" %4.3f",values[knn[i]]);
  printf("\n");
  printf(" knn:");
  for (int i=0;i<k;i++)
    printf(" %3i",knn[i]);
  printf("\n");
  }
#endif

  if (sumw < LOW_WEIGHT)
  {
    //if (pt==151)
    //printf(" low weights \n");
    nlw++;

    delete [] knn;
    delete [] knnd;
    delete [] knnw;

    return 1.2;
  }

  for (int i=0;i<k;i++)
    val += knnw[i] * values[knn[i]];


  delete [] knn;
  delete [] knnd;
  delete [] knnw;

  return val;
}


int KNNR::find_knn(int pt, int k, int* knn, double* knnd)
{
  double* X1 = new double[fsize];
  for (int i=0;i<fsize;i++)
    X1[i] = X[pt*fsize+i];

  ptskip = pt;
  int nfound = find_knn(X1,k,knn,knnd);
  ptskip = -1;

  delete [] X1;

  return nfound;
}

int KNNR::check_unique(int pt)
{
  if (udata==NULL)
    return 1;
  else if (udata[pt]==0)
    return 0;

  return 1;
}

int KNNR::find_knn(double* X1, int k, int* knn, double* knnd)
{
  int nfound = 0;
  if (k>knn_N)
    return 0;

  int* close_n = knn;
  double* close_d = knnd;
  for (int i=0;i<k;i++) close_n[i] = -1;
  for (int i=0;i<k;i++) close_d[i] = 10000.;

  for (int i=0;i<knn_N;i++)
  if (i!=ptskip && check_unique(i))
  {
    double dist1 = get_distance(X1,&X[i*fsize]);
    for (int j=0;j<k;j++)
    {
     //important: find k nearest neighbors, priority to high values
      double valcomp = 0.;
      if (close_n[j]!=-1) valcomp = values[close_n[j]];

      int moveup = 0;
//      if (dist1<close_d[j]-CZERO) moveup = 1;
      if (dist1<close_d[j]) moveup = 1;
//was on, was using CZERO (bad)
      if (values[i]>valcomp && close_val(dist1,close_d[j],DZERO)) moveup = 1;
      if (moveup)
      {
        //printf(" %4.3f is closer than %4.3f \n",dist1,close_d[j]);
        //printf(" close_d: %4.3f %4.3f %4.3f \n",close_d[0],close_d[1],close_d[2]);
        for (int l=k-1;l>j;l--)
        {
          close_d[l] = close_d[l-1];
          close_n[l] = close_n[l-1];
        }
        close_d[j] = dist1;
        close_n[j] = i;
  
        //printf(" close_d: %4.3f %4.3f %4.3f \n",close_d[0],close_d[1],close_d[2]);
  
        break;
      }
    } //loop j over current close set
  } //loop i over npts
  
  //printf(" final close_d: %4.3f %4.3f %4.3f \n",close_d[0],close_d[1],close_d[2]);
  //printf(" final close_n: %3i %3i %3i \n",close_n[0],close_n[1],close_n[2]);

  for (int i=0;i<k;i++)
  if (close_n[i]>-1)
    nfound++;


  return nfound;
}


void KNNR::get_distances(int npts1, double* X0)
{
  if (X0==NULL) return;

  //printf(" in get_distances npts1: %i \n",npts1); fflush(stdout);
  if (distances!=NULL) 
    delete [] distances;
  distances = new double[npts1*npts1];

  for (int n=0;n<npts1;n++)
  for (int m=0;m<n;m++)
  {
    distances[n*npts1+m] = get_distance(&X0[n*fsize],&X0[m*fsize]);
  }
  for (int n=0;n<npts1;n++)
  for (int m=n+1;m<npts1;m++)
    distances[n*npts1+m] = distances[m*npts1+n];


printf(" in get_distances, quiet: %i \n",quiet);
  if (quiet)
    return;

#if PRINT_DISTS
  if (0)
  if (!quiet)
  for (int n=0;n<npts1;n++)
  {
    for (int m=0;m<npts1;m++)
      printf(" %5.4f",distances[n*npts1+m]);
    printf("\n");
  }
#endif

#if 1
  int* ids1 = new int[npts1];
  for (int i=0;i<npts1;i++)
    ids1[i] = ids[i];
  for (int i=0;i<npts1;i++)
  //if (!udata[i])
  if (!check_unique(i))
    ids1[i] += 1000000;

 //gephi data
//  double ALPHA_SAVE = KNNR_ALPHA;
  double ALPHA_SAVE = 0.15;
  double SAVE_THRESH = 0.5;
  int nf = 0;
  for (int i=0;i<npts1;i++)
  for (int j=0;j<i;j++)
  if (exp(-distances[i*npts1+j]/ALPHA_SAVE)<SAVE_THRESH)
    nf++;
//  if (nf<5000)
    SAVE_THRESH /= 10.;
//  save_gephi(npts1,ids1,distances,values_print,ALPHA_SAVE,SAVE_THRESH);
  save_gephi_3(npts1,ids1,distances,values,values_print,ALPHA_SAVE,SAVE_THRESH);

  delete [] ids1;
#endif

  return;
}

double KNNR::get_distance(double* X1, double* X2)
{
  double d = 0.;

  double diff;
  for (int i=0;i<fsize;i++)
  if (active[i])
  {
    diff = X1[i] - X2[i];
    d += diff*diff;
  }
  d = sqrt(d);

  return d;
}

void KNNR::load_values_print(int npts1, double* y1)
{
  if (values_print!=NULL)
    delete [] values_print;
  values_print = new double[npts];

  for (int i=0;i<npts;i++)
    values_print[i] = y1[i];

  return;
}

void KNNR::load_values(int npts1, int fsize1, double* X1, double* y1)
{
  if (!quiet)
    printf(" in KNNR load_values for %i pts, fsize: %i \n",npts1,fsize1);

  npts = npts1;
  knn_N = npts;
  fsize = fsize1;

  if (values!=NULL)
    delete [] values;
  values = new double[npts];

  load_values_print(npts1,y1);
  for (int i=0;i<npts;i++)
    values[i] = y1[i];

  if (X!=NULL)
    delete [] X;
  X = new double[npts*fsize];

  for (int i=0;i<npts*fsize;i++)
    X[i] = X1[i];

  determine_active();

#if 0
  printf("\n printing X,y: \n");
  for (int i=0;i<npts;i++)
  {
    for (int j=0;j<fsize;j++)
      printf(" %2.1f",X[i*fsize+j]);
    printf("  %4.3f   id: %3i/%3i \n",values[i],i,ids[i]);
  }
  printf("\n");
#endif


  return;
}

void KNNR::determine_active()
{
//  if (!quiet)
//    printf("  in determine_active: fsize: %i \n",fsize);

  if (active!=NULL)
    delete [] active;
  active = new double[fsize];

  for (int i=0;i<fsize;i++)
    active[i] = 0.;

  double acv = 0.;
  for (int i=0;i<fsize;i++)
  {
    acv = X[0*fsize+i];
    for (int j=1;j<npts;j++)
    {
      if (!close_val(X[j*fsize+i],acv,0.0001))
      {
        active[i] = 1.;
        break;
      }
    }
  }

  return;
}

void KNNR::freemem()
{
  if (values!=NULL)
    delete [] values;
  if (distances!=NULL)
    delete [] distances;
  if (X!=NULL)
    delete [] X;
}

void KNNR::reset_active()
{
  if (active!=NULL)
    delete [] active;
  active = NULL;
}

void KNNR::init()
{
  quiet = 0;

  pthresh = 0.5;
  LOW_WEIGHT = 0.1;

  nlw = 0;
  ptskip = -1;
  npts = 0;
  knn_N = 0;
  fsize = 0;
  values = NULL;
  values_print = NULL;
  X = NULL;

  distances = NULL;
  active = NULL;
  ids = NULL;
  udata = NULL;

  errlist = NULL;
  knnlist = NULL;

  return;
}



