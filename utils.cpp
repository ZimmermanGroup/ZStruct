// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "utils.h"
#include "omp.h"
using namespace std;
//#include <Accelerate/Accelerate.h>
//#include <vecLib/clapack.h>
//#include "/export/apps/Intel/Compiler/11.1/075/mkl/include/mkl.h"
#include <mkl.h>

const double BOHRtoANG=0.52917720859;
const double ANGtoBOHR =1.0000000/BOHRtoANG;



//sorts v1, arranging v2 with it (data in pairs)
void sort_n2p(double* v1, double* v2)
{
  double v1a1 = v1[0];
  if (v1a1<v1[2])
  {
    double v1a2 = v1[1];
    double v2a1 = v2[0];
    double v2a2 = v2[1];
    v1[0] = v1[2];
    v1[1] = v1[3];
    v1[2] = v1a1;
    v1[3] = v1a2;
    v2[0] = v2[2];
    v2[1] = v2[3];
    v2[2] = v2a1;
    v2[3] = v2a2;
  }
  else if (v1a1==v1[2] && v1[1]<v1[3])
  {
    double v1a2 = v1[1];
    double v2a1 = v2[0];
    double v2a2 = v2[1];
    v1[0] = v1[2];
    v1[1] = v1[3];
    v1[2] = v1a1;
    v1[3] = v1a2;
    v2[0] = v2[2];
    v2[1] = v2[3];
    v2[2] = v2a1;
    v2[3] = v2a2;
  }

  return;
}

void sort_n3p(double* v1, double* v2)
{
  sort_n2p(&v1[2],&v2[2]);
  sort_n2p(v1,v2);
  sort_n2p(&v1[2],&v2[2]);

  return;
}

void sort_n4p(double* v1, double* v2)
{
  sort_n3p(&v1[2],&v2[2]);
  sort_n3p(v1,v2);
  sort_n3p(&v1[2],&v2[2]);

  return;
}

void sort_n5p(double* v1, double* v2)
{
  sort_n4p(&v1[2],&v2[2]);
  sort_n4p(v1,v2);
  sort_n4p(&v1[2],&v2[2]);

  return;
}

//decreasing order
void sort_2(double* v1)
{
  double v1a = v1[0];
  if (v1a<v1[1])
  {
    v1[0] = v1[1];
    v1[1] = v1a;
  }

  return;
}

void sort_3(double* v1)
{
  sort_2(&v1[1]);
  sort_2(v1);
  sort_2(&v1[1]);

  return;
}

void sort_4(double* v1)
{
  sort_3(&v1[1]);
  sort_3(v1);
  sort_3(&v1[1]);
  return;
}

void sort_5(double* v1)
{
  sort_4(&v1[1]);
  sort_4(v1);
  sort_4(&v1[1]);
  return;
}

//sorts v1, arranging v2 with it
void sort_2p(double* v1, double* v2)
{
  double v1a = v1[0];
  if (v1a<v1[1])
  {
    double v2a = v2[0];
    v1[0] = v1[1];
    v1[1] = v1a;
    v2[0] = v2[1];
    v2[1] = v2a;
  }
  else if (v1a==v1[1] && v2[0]<v2[1])
  {
    double v2a = v2[0];
    v1[0] = v1[1];
    v1[1] = v1a;
    v2[0] = v2[1];
    v2[1] = v2a;
  }

  return;
}

void sort_3p(double* v1, double* v2)
{
  sort_2p(&v1[1],&v2[1]);
  sort_2p(v1,v2);
  sort_2p(&v1[1],&v2[1]);

  return;
}

void sort_4p(double* v1, double* v2)
{
  sort_3p(&v1[1],&v2[1]);
  sort_3p(v1,v2);
  sort_3p(&v1[1],&v2[1]);
  return;
}

void sort_5p(double* v1, double* v2)
{
  sort_4p(&v1[1],&v2[1]);
  sort_4p(v1,v2);
  sort_4p(&v1[1],&v2[1]);
  return;
}

double sigmoid(double A, double	A0, double alpha)
{
  double term = (A-A0)/alpha;
  return 1/(1+exp(term));
}

double norm(double* x, int size)
{
  double val = 0.;
  for (int i=0;i<size;i++)
    val += x[i]*x[i];
  val = sqrt(val);

  return val;
}

int check_array(int size, double* A)
{
  int bad = 0;
  for (int i=0;i<size;i++)
  if (A[i]!=A[i])
  {
    bad = 1;
    break;
  }

  return bad;
}

double read_chk_file(string filename, int natoms, double* xyz1)
{
  ifstream output(filename.c_str(),ios::in);
  if (!output) { printf(" error opening chk file: %s \n",filename.c_str()); return -999.; }

  double energy = 999.;
  string line;
  vector<string> tok_line;

  int count = 0;
  while(!output.eof())
  {
    getline(output,line);
    //cout << " RR " << line << endl;
    if (line.find("Current cartesian coordinates")!=string::npos)
    {
      int nf = 0;
      while (nf<3*natoms)
      {
        getline(output,line);
        tok_line = StringTools::tokenize(line, " \t");
        int lsize = tok_line.size();
        for (int j=0;j<lsize;j++)
          xyz1[nf++] = atof(tok_line[j].c_str())*BOHRtoANG;
      }
      count++;
    } 
    if (line.find("SCF Energy")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      energy = atof(tok_line[3].c_str());
      count++;
    }

    if (count>1) break;
  }

  output.close();

#if 1
  for (int i=0;i<natoms;i++)
    printf(" %6.5f %6.5f %6.5f \n",xyz1[3*i+0],xyz1[3*i+1],xyz1[3*i+2]);
#endif

  return energy;
}

double read_xyz_file(string filename, int natoms, double* xyz1)
{
  ifstream output(filename.c_str(),ios::in);
  if (!output) { printf(" error opening xyz file: %s \n",filename.c_str()); return -999.; }

  string line;
  vector<string> tok_line;

  getline(output,line);
  int natoms1 = atoi(line.c_str());
  if (natoms!=natoms1)
  {
    printf(" natoms mismatch, exiting!  (%i/%i) \n",natoms,natoms1);
    exit(1);
  }

  getline(output,line);
  double energy = atof(line.c_str());

  int count = 0;
  while(!output.eof())
  {
    getline(output,line);
    tok_line = StringTools::tokenize(line, " \t[]");
 
    //cout << " RR " << line << endl;

    if (tok_line.size()<3)
    {
      printf("\n ERROR: bad xyz file \n");
      exit(1);
    }

    if (tok_line.size()==3)
    {
      xyz1[3*count+0] = atof(tok_line[0].c_str());
      xyz1[3*count+1] = atof(tok_line[1].c_str());
      xyz1[3*count+2] = atof(tok_line[2].c_str());
    }
    else if (tok_line.size()==4)
    {
      xyz1[3*count+0] = atof(tok_line[1].c_str());
      xyz1[3*count+1] = atof(tok_line[2].c_str());
      xyz1[3*count+2] = atof(tok_line[3].c_str());
    }

    count++;
    if (count>=natoms) break;
  }

  output.close();

#if 1
  for (int i=0;i<natoms;i++)
    printf(" %6.5f %6.5f %6.5f \n",xyz1[3*i+0],xyz1[3*i+1],xyz1[3*i+2]);
#endif

  return energy;
}

void get_rotation_matrix(double** rotMat, double* thetas)
{
  double x=thetas[0]; double y=thetas[1]; double z=thetas[2];
  rotMat[0][0] = cos(y)*cos(z);
  rotMat[0][1] = -cos(y)*sin(z);
  rotMat[0][2] = sin(y);
  rotMat[1][0] = sin(x)*sin(y)*cos(z)+cos(x)*sin(z);
  rotMat[1][1] = -sin(x)*sin(y)*sin(z)+cos(x)*cos(z);
  rotMat[1][2] = -sin(x)*cos(y);
  rotMat[2][0] = -cos(x)*sin(y)*cos(z)+sin(x)*sin(z);
  rotMat[2][1] = cos(x)*sin(y)*sin(z)+sin(x)*cos(z);
  rotMat[2][2] = cos(x)*cos(y);
}

//A = B*C
int matmult_sq(double* A, double* B, double* C, int N)
{
  for (int i=0;i<N*N;i++) A[i] = 0.;
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    A[i*N+j] += B[i*N+k] * C[k*N+j];

  return 0;
}

int fact(int n)
{
  int a = n;
  for (int i=1;i<n;i++)
    a *= n-i;
  return a;
}

int randomi(int a){

  timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec*t1.tv_sec);

  int randn = rand()%a;
  return randn;
}

double randomf(double a)
{

  timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec*t1.tv_sec);

  double r = (rand()%100000000)/100000000.; //8 digits of rand
  //printf(" r: %4.3f \n",r);
  double randn = r*a;
  return randn;
}

int sign(double x){

  if (x>=0) return 1;
  return -1;

}

void cross(double* m, double* r1, double* r2){

  m[0] = r1[1]*r2[2]-r2[1]*r1[2];
  m[1] = -r1[0]*r2[2]+r2[0]*r1[2];
  m[2] = r1[0]*r2[1]-r2[0]*r1[1];

  return;
}

void trans(double* Bt, double* B, int m, int n) {

  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
    Bt[i*n+j] = B[j*m+i];

  return;
}

int close_val(double x1, double x2, double diff)
{
  int close = 0;
  if (fabs(x1-x2)<diff)
    close = 1;
  return close;
}


int Invert(double* A, int m){

  int error = 0;

#if 0
  printf(" in Invert call, m: %i \n",m);
  for (int i=0;i<m;i++)
  {
    for (int j=0;j<m;j++)
      printf(" %1.3f",A[m*i+j]);
    printf("\n");
  }
  fflush(stdout);
#endif

  int LenWork = 4*m;
  double* Work = new double[LenWork];

  int Info = 0;

  //printf(" LenWork: %i \n",LenWork);

  int* IPiv = new int[m];

//printf("\n");
//  dgesdd_((char*)"A", &m, &n, A, &LDA, S, U, &m, Vt, &n, Work, &LenWork, IWork, &Info);

  dgetrf_(&m,&m,A,&m,IPiv,&Info);
  if (Info!=0)
    printf(" after dgetrf, Info error is: %i \n",Info);

  dgetri_(&m,A,&m,IPiv,Work,&LenWork,&Info);
  if (Info!=0)
  {
    printf(" after invert, Info error is: %i \n",Info);
    printf(" A-1: \n");
    for (int i=0;i<m;i++)
    {
      for (int j=0;j<m;j++)
        printf(" %4.3f",A[i*m+j]);
      printf("\n");
    }
    error = 1;
  }


  delete [] IPiv;
  delete [] Work;


  return error;
}


int Diagonalize(double* A, double* eigen, int size){

#ifdef _OPENMP
//  mkl_set_num_threads(1);
//  int nthreads = omp_get_num_threads();
//  omp_set_num_threads(1);
#endif

#define DSYEVX 1
 // printf(" in diagonalize call, size: %i \n",size);
 // printf(" in diagonalize: mkl_threads: %i \n",mkl_get_max_threads());

  int N = size;
  int LDA = size;
  double* EVal = eigen;

//borrowed from qchem liblas/diagon.C


    char JobZ = 'V', Range = 'A', UpLo = 'U';
    int IL = 1, IU = N;
    double AbsTol = 0.0, VL = 1.0, VU = -1.0;
    int NEValFound;

    double* EVec = new double[LDA*N];

    // Give dsyevx more work space than the minimum 8*N to improve performance
#if DSYEVX
    int LenWork = 32*N; //8*N min for dsyevx
#else
    int LenWork = 1+6*N+2*N*N; //1+6*N+2*N*N min for dsyevd
#endif
    double* Work = new double[LenWork];

#if DSYEVX
    int LenIWork = 5*N; //5*N for dsyevx
#else
    int LenIWork = 10*N; //3+5*N min for dsyevd
#endif
    int* IWork = new int[LenIWork];
    int* IFail = new int[N];

    int Info = 0;

#if DSYEVX
    dsyevx_(&JobZ, &Range, &UpLo, &N, A, &LDA, &VL, &VU, &IL, &IU, &AbsTol,
           &NEValFound, EVal, EVec, &LDA, Work, &LenWork, IWork, IFail, &Info);
#else
    dsyevd_(&JobZ, &UpLo, &N, A, &LDA, EVal, Work, &LenWork, IWork, &LenIWork, &Info);
#endif

#if 0
    if (Info != 0 && KillJob) {
      printf(" Info = %d\n",Info);
      QCrash("Call to dsyevx failed in Diagonalize");
    }
#endif

  int n_nonzero = 0;
  for (int i=0;i<size;i++)
  {
    //printf(" eigenvalue %i: %1.5f \n",i,eigen[i]);
    if (abs(eigen[i])>0.0001) n_nonzero++;
  }
  //printf(" found %i independent vectors \n",n_nonzero);

#if DSYEVX
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    A[i*size+j] = EVec[i*size+j];
#endif

#ifdef _OPENMP
//  omp_set_num_threads(nthreads);
#endif

    delete [] EVec;
    delete [] Work;
    delete [] IWork;
    delete [] IFail;

  return 0;

}

int Diagonalize(double* A, double* eigenvecs, double* eigen, int size){

  printf(" in diagonalize call, size: %i \n",size);
  

  int N = size;
  int LDA = size;
  double* EVal = eigen;

//borrowed from qchem liblas/diagon.C


    char JobZ = 'V', Range = 'A', UpLo = 'U';
    int IL = 1, IU = N;
    double AbsTol = 0.0, VL = 1.0, VU = -1.0;
    int NEValFound;

    double* EVec = new double[LDA*N];

    // Give dsyevx more work space than the minimum 8*N to improve performance
    int LenWork = 32*N;
    double* Work = new double[LenWork];

    int LenIWork = 5*N;
    int* IWork = new int[LenIWork];
    int* IFail = new int[N];

    int Info = 0;

    dsyevx_(&JobZ, &Range, &UpLo, &N, A, &LDA, &VL, &VU, &IL, &IU, &AbsTol,
           &NEValFound, EVal, EVec, &LDA, Work, &LenWork, IWork, IFail, &Info);

#if 0
    if (Info != 0 && KillJob) {
      printf(" Info = %d\n",Info);
      QCrash("Call to dsyevx failed in Diagonalize");
    }
#endif

  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    eigenvecs[i*size+j] = EVec[i*size+j];

    delete [] EVec;
    delete [] Work;
    delete [] IWork;
    delete [] IFail;

  return 0;

}

void order_array(int* olist, double* Ao, double* A, int na)
{
#if 0
  printf(" before: \n");
  for (int i=0;i<na;i++)
    printf(" %i: %4.3f\n",i,A[i]);
#endif
  for (int i=0;i<na;i++)
  {
    olist[i] = -1;
    Ao[i] = 100000;
  }
  for (int i=0;i<na;i++)
  {
    int nmax = 0;
    double max = -100000;
    for (int j=0;j<na;j++)
    if (A[j]>max)
    {
      int found = 0;
      for (int k=0;k<na;k++)
      if (olist[k]==j)
        found = 1;
      if (!found)
      {
        nmax = j;
        max = A[j];
      }
    } //loop j over search for i'th entry
    olist[i] = nmax;
    Ao[i] = max;
  } //loop i 

#if 0
  printf(" ordered: \n");
  for (int i=0;i<na;i++)
    printf(" %i: %4.3f\n",i,Ao[i]);
  printf(" olist:");
  for (int i=0;i<na;i++)
    printf(" %i",olist[i]);
  printf("\n");
#endif

  return;
}

double getRa(int* anumbers, int a1)
{
  double value;
  if	  (anumbers[a1]==1) value = 1.3;
  else if (anumbers[a1]==3) value = 2.65; //PT
  else if (anumbers[a1]==4) value = 2.0; //PT
  else if (anumbers[a1]==5) value = 1.75;
  else if (anumbers[a1]==6) value = 1.65;
  else if (anumbers[a1]==7) value = 1.65;
  else if (anumbers[a1]==8) value = 1.65;
  else if (anumbers[a1]==9) value = 1.6;
  else if (anumbers[a1]==11) value = 3.3; //PT
  else if (anumbers[a1]==12) value = 3.1;
  else if (anumbers[a1]==13) value = 2.6;
  else if (anumbers[a1]==14) value = 2.6;
  else if (anumbers[a1]==15) value = 2.5;
  else if (anumbers[a1]==16) value = 2.45;
  else if (anumbers[a1]==17) value = 2.1;
  else if (anumbers[a1]==19) value = 4.0; //PT
  else if (anumbers[a1]==20) value = 3.5; //PT
  else if (anumbers[a1]==21) value = 4.5;
  else if (anumbers[a1]==22) value = 4.2;
  else if (anumbers[a1]==23) value = 4.0;
  else if (anumbers[a1]==24) value = 3.5;
  else if (anumbers[a1]==25) value = 3.4;
  else if (anumbers[a1]==26) value = 3.3;
  else if (anumbers[a1]==27) value = 3.0;
  else if (anumbers[a1]==28) value = 3.0;
  else if (anumbers[a1]==29) value = 3.0;
  else if (anumbers[a1]==30) value = 3.0;
  else if (anumbers[a1]==35) value = 2.7;
  else if (anumbers[a1]==40) value = 3.35;
  else if (anumbers[a1]==44) value = 3.2;
  else if (anumbers[a1]==45) value = 3.15;
  else if (anumbers[a1]==46) value = 3.15;
  else if (anumbers[a1]==47) value = 3.25;
  else if (anumbers[a1]==53) value = 2.8; //iodine
  else if (anumbers[a1]==73) value = 3.3;
  else if (anumbers[a1]==74) value = 3.3;
  else if (anumbers[a1]==75) value = 3.3;
  else if (anumbers[a1]==76) value = 3.3;
  else if (anumbers[a1]==77) value = 3.35;
  else if (anumbers[a1]==78) value = 3.35;
  else if (anumbers[a1]==79) value = 3.45;
  else if (anumbers[a1]==0) value = 2.5;
  else
  {
    printf(" Need to add atomic number %i to getR! \n",anumbers[a1]);
    exit(1);
  }

  return value;
}
