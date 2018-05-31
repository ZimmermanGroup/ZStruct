// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include <iostream>
#include <fstream>
#include <stdio.h>

#include "zstruct.h"
#include "rxndb.h"
#include "knnr.h"
#include "lsq.h"
#include "rxnftr.h"
#include "icoord.h"
#include "align.h"
#include "utils.h"
double read_temperature();

using namespace std;


int main(int argc, char* argv[])
{

  string inpfile;
  string xyzfile;
  double emax = 250.;
  double emin = -1000.;
  int nsteps = 1;
  xyzfile="test.xyz";
  switch (argc){
  case 1:
//    printf(" case 1 \n");
    break;
  case 2:
//    printf(" case 2 \n");
    nsteps = atoi(argv[1]);
    break;
  case 3:
//    printf(" case 3 \n");
    emax = atof(argv[1]);
    nsteps = atoi(argv[2]);
    break;
  case 4:
//    printf(" case 4 \n");
    emax = atof(argv[1]);
    emin = atof(argv[2]);
    nsteps = atoi(argv[3]);
    break;
  default:
    inpfile="xyzfile";
    break;
//    return -1;
  }

  int done=0;




  if (nsteps==-1)
  {
    printf("  restarting regression from file! \n");
    RXNDB rxndb1;
    double pthresh = 0.4;

    KNNR knnr1;
    knnr1.init();
    knnr1.quiet = 0;
    knnr1.pthresh = pthresh;

    int N = 10000;
    int dim = 40;
    double* X = new double[N*dim];
    double* Ea = new double[N];
    double* Erxn = new double[N];
    int* ids = new int[N];
    N = rxndb1.read_xy_data(X,Ea,Erxn,ids,"xydata");
    knnr1.ids = ids;

    double* Pr = new double[N];
    double* PrT2 = new double[N];  
    double temperature = read_temperature();
    double kT = temperature*(0.6/300.); //need precision

#if 1
    double LOW_WEIGHT = 0.000;
    double LOW_WEIGHT_INTERVAL = 0.02;
    double TRAINSHIFT = 0.;
    double TRAIN_INTERVAL = 5.0;
#else
    double LOW_WEIGHT = 0.00;
    double LOW_WEIGHT_INTERVAL = 0.02;
    double TRAINSHIFT = 10.;
    double TRAIN_INTERVAL = 5.0;
#endif

    double Eref = 60*kT;
    double alpha_s = 10.*kT; //was 5kT

#if 0
    for (int i=0;i<N;i++)
      Pr[i] = sigmoid(Ea[i],Eref,alpha_s);
#else
    printf("  using EA instead of Pr! \n");
    for (int i=0;i<N;i++) Pr[i] = Ea[i];
#endif

    printf(" T: %4.2f K  kT: %3.2f kcal/mol  Eref: %3.2f kcal/mol  sigmoid_alpha: %3.2f  N: %2i \n",temperature,kT,Eref,alpha_s,N);
    for (int n=0;n<1;n++)
    {
      printf("\n\n  Doing TRAINSHIFT: %4.1f LOW_WEIGHT: %5.2f \n\n",TRAINSHIFT,LOW_WEIGHT);

      knnr1.LOW_WEIGHT = LOW_WEIGHT;
      for (int i=0;i<N;i++)
        PrT2[i] = sigmoid(Ea[i],Eref+TRAINSHIFT,alpha_s);
      for (int i=0;i<N;i++) PrT2[i] = Ea[i];

      knnr1.load_values(N,dim,X,PrT2);
      knnr1.load_values_print(N,Pr);

      int KNNR_K = 4;
      double error = knnr1.test_points(KNNR_K);

      string scpstr = "scp";
      printf("\n printing false negatives \n");
      int nfalseneg = 0;
      for (int i=0;i<N;i++)
      {
        double val = knnr1.errlist[i] + PrT2[i];
        if (val<Pr[i]/2.0 && Pr[i]>pthresh)
        {
          printf("  pt: %4i id: %4i val: %3.2f actual: %3.2f Ea: %5.1f \n",i,ids[i],val,Pr[i],Ea[i]);
          for (int j=0;j<dim;j++)
            printf(" %4.1f",X[i*dim+j]*10.);
          printf("\n");
          for (int l=0;l<KNNR_K;l++)
          {
            int index = knnr1.knnlist[i*KNNR_K+l];
   	    for (int j=0;j<dim;j++)
              printf(" %4.1f",X[index*dim+j]*10.);
            printf("   %4.3f %4.3f  id: %5i/%5i",Pr[index],PrT2[index],index,ids[index]);
            printf("  dist: %4.3f \n",knnr1.get_distance(&X[i*dim],&X[index*dim]));
          }
          string nstr = StringTools::int2str(ids[i],4,"0");
          scpstr += " stringfile.xyz"+nstr;
          nfalseneg++;
        }
      }
      printf("  found %3i false negatives \n",nfalseneg);
      printf(" %s $I \n",scpstr.c_str());

      TRAINSHIFT += TRAIN_INTERVAL;
//      LOW_WEIGHT += LOW_WEIGHT_INTERVAL;
    } //loop n over bias

    delete [] Pr;
    delete [] PrT2;

    exit(1);
  }







  if (nsteps>10) 
  {
    printf(" nsteps>10 is (maybe) a bad idea \n");
    nsteps = 10;
  }
  printf(" (main) emax is: %1.4f \n",emax);
  printf(" (main) emin is: %1.4f \n",emin);
  printf(" (main) nsteps is: %i \n",nsteps);

  printf("\n\n\n");


#if 1
  ZStruct zmain;
//  zmain.init(xyzfile,xyzlist);
  zmain.go_zstruct(1);
#endif

#if 0
  ICoord ic1,ic2;
  ic1.init("react1.xyz");
  ic2.init("react2.xyz");
  printf("\n");

  Align align1;
  align1.init(ic1.natoms,ic1.anames,ic1.anumbers,ic1.coords,ic2.natoms,ic2.anames,ic2.anumbers,ic2.coords);

  int nadd = 2;
  int* add = new int[4];
  add[0] = 0; add[1] = 6; add[2] = 1; add[3] = 7;
  align1.add_align(nadd,add);
  align1.print_xyz();

#endif  

#if 0
  LSQ lsq1;
  lsq1.init();
  lsq1.go_test();
#endif

#if 0
  int natoms = 3;
  RXNFTR rxnftr1;
  rxnftr1.init(natoms);
  rxnftr1.set_element(0,1);
  rxnftr1.set_abs(0,-1,-1,1);
  rxnftr1.set_element(1,6);
  rxnftr1.set_coordn(1,3);
  rxnftr1.print_atoms();
#endif

  return 0;
}



double read_temperature()
{
  string filename = "TEMPERATURE";

  ifstream file;
  file.open(filename.c_str());
  if (!file)
  {
    printf(" couldn't find %s file \n",filename.c_str());
    return -1.;
  }

  string line;
  getline(file, line);
  int length=StringTools::cleanstring(line);
  double T = atof(line.c_str());

  file.close();

//  printf(" found T: %4.2f \n",T);
    
  return T;
}

