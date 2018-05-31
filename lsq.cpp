// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "lsq.h"
using namespace std;


//do ridge regression?


void LSQ::init()
{
  params = NULL;
  quiet = 0;

  return;
}

void LSQ::go_test()
{
  printf("\n in LSQ go_test \n");

  double RANDNOISE = 0.1;

  N = 10;
  double* y = new double[N];
  for (int i=0;i<N;i++)
    y[i] = i/2.0 + 0.0 + RANDNOISE/2.0-randomf(RANDNOISE);

  p = 2;
  double* x = new double[p*N];
  for (int i=0;i<N;i++)
    x[p*i+0] = i;
  if (p==2)
  for (int i=0;i<N;i++)
    x[p*i+1] = 0.;

  printf(" x: \n");
  for (int i=0;i<N;i++)
  {
    for (int j=0;j<p;j++)
      printf(" %4.3f",x[p*i+j]);
    printf("\n");
  }

  printf(" y: \n");
  for (int i=0;i<N;i++)
    printf(" %4.3f\n",y[i]);
  printf("\n");

  double error = do_least_squares(x,p,y,N);

  printf("   total unsigned error: %4.3f \n",error);
  printf(" average unsigned error: %4.3f \n",error/N);

  delete [] x;
  delete [] y;

  return;
}


double LSQ::do_least_squares(double* x0, int p0, double* y, int N0)
{
  p = p0;
  N = N0;

  if (!quiet)
    printf(" creating regression (p: %i N: %i) \n",p,N);

  double* x = x0;

  double* xT = new double[p*N];
  trans(xT,x,p,N);

#if 0
  printf(" xT: \n");
  for (int i=0;i<p;i++)
  {
    for (int j=0;j<N;j++)
      printf(" %4.3f",xT[i*N+j]);
    printf("\n");
  }
#endif
 
  double* xTx = new double[p*p];
  for (int i=0;i<p*p;i++) xTx[i] = 0.;
  for (int i=0;i<p;i++)
  for (int j=0;j<p;j++)
  for (int k=0;k<N;k++)
    xTx[i*p+j] += xT[i*N+k] * x[k*p+j];

#if 0
  printf(" xTx: \n");
  for (int i=0;i<p;i++)
  {
    for (int j=0;j<p;j++)
      printf(" %4.3f",xTx[i*p+j]);
    printf("\n");
  }
#endif

  //adjust zero diagonals so invert works
  for (int i=0;i<p;i++)
  if (fabs(xTx[i*p+i])<0.000001)
    xTx[i*p+i] = 10000.;

#if 0
  //debug
  for (int i=0;i<p*p;i++)
    xTx[i] = 0.;
  for (int i=0;i<p;i++)
    xTx[i*p+i] = 1.;
#endif

#if 1
  //debug
  double* tmp = new double[p*p];
  for (int i=0;i<p*p;i++) tmp[i] = xTx[i];
  double* eig = new double[p];
  for (int i=0;i<p;i++) eig[i] = 0.;
  Diagonalize(tmp,eig,p);
  if (!quiet)
  {
    printf("\n eigenvalues:");
    for (int i=0;i<p;i++)
      printf(" %8.6f",eig[i]);
    printf("\n");
  }


  int zeroeig = 0;
  for (int i=0;i<p;i++)
  if (fabs(eig[i])<0.0000001)
    zeroeig++;
  if (zeroeig)
  {
    if (!quiet)
    {
      printf(" found zero eigenvalue(s) (%i), cannot invert \n",zeroeig);
      printf(" p: %i N: %i \n",p,N);
    }
    delete [] xT;
    delete [] xTx;
    delete [] tmp;
    delete [] eig;
    return 99999.;
  }

  delete [] tmp;
  delete [] eig;
#endif

  int error1 = 0;
  if (p>1) error1 = Invert(xTx,p);
  else xTx[0] = 1./xTx[0];
  if (error1)
  {
    printf("\n failed to invert xTx \n");
    return 99999.;
  }

#if 0
  printf(" xTxi: \n");
  for (int i=0;i<p;i++)
  {
    for (int j=0;j<p;j++)
      printf(" %4.3f",xTx[i*p+j]);
    printf("\n");
  }
  fflush(stdout);
#endif

  double* A = new double[p*N];
  for (int i=0;i<p*N;i++) A[i] = 0.;
  for (int i=0;i<p;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<p;k++)
    A[i*N+j] += xTx[i*p+k] * xT[k*N+j];

#if 0
  printf(" A: \n");
  for (int i=0;i<p;i++)
  {
    for (int j=0;j<N;j++)
      printf(" %4.3f",A[i*N+j]);
    printf("\n");
  }
#endif

  if (params!=NULL) delete [] params;
  params = new double[p];
  for (int i=0;i<p;i++) params[i] = 0.;
  for (int i=0;i<p;i++)
  for (int k=0;k<N;k++)
    params[i] += A[i*N+k] * y[k];

  if (!quiet)
  {
    printf("\n params:");
    for (int i=0;i<p;i++)
      printf(" %6.4f",params[i]);
    printf("\n");
    fflush(stdout);
  }

  if (!quiet)
  {
    int non_contribute = 0;
    for (int i=0;i<p;i++)
    if (close_val(params[i],0.,0.00001))
      non_contribute++;
    printf(" %i parameters doing nothing \n\n",non_contribute);
  }

  double error = 0.;
  double err,val;
  for (int i=0;i<N;i++)
  {
    val = estimate_point(&x[p*i]);
#if 1
    err = val-y[i];
    error += err*err*(1+y[i]); //weighted error
#elif 0
    err = val-y[i];
    error += err*err;
#else
    error += fabs(val-y[i]);
#endif
    if (!quiet)
      printf(" actual pred: %6.3f %6.3f \n",y[i],val);
  }
  error = sqrt(error / N);

  delete [] A;
  delete [] xTx;
  delete [] xT;

  return error;
}


double LSQ::estimate_point(double* x)
{
  double val = 0.;

  for (int j=0;j<p;j++)
    val += x[j] * params[j];

  return val;
}
