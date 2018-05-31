// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "icoord.h"
#include "utils.h"

//for CG optimizer
#define SCALESD 25
#define SCALE1 1.4
#define SCALEA 3.0 
#define OPTTHRESH 0.01
#define MAXAD 0.15

#define OPTSTEPS 500
#define SHADOWSTEPS 40
 
int ICoord::opt(){

  string xyzfile = "scratch/xyzfile.xyz";
  return opt(xyzfile);
}

int ICoord::opt(string xyzfile_string){

//  printf("  \n"); 
  dxm1 = new double[3*natoms];

  mm_grad();
  //print_grad();

  ofstream xyzfile;
//  string xyzfile_string = "xyzfile.txt";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

  for (int n=0;n<OPTSTEPS;n++)
  {
//    printf(" Opt step: %i ",n+1);

    xyzfile << " " << natoms << endl;
    xyzfile << " " << gradrms << endl;
    for (int i=0;i<natoms;i++) 
    {
      xyzfile << "  " << anames[i];
      xyzfile << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2];
      xyzfile << endl;
    }

    if (n==30) 
      pgradrms = 10000; // reset CG optimizer
    if (n%100==0) 
      pgradrms = 10000; // reset CG optimizer
    if (n==0) update_xyz_sd();
    else update_xyz_cg();
    mm_grad();
//    print_grad();
//    printf(" gradient RMS: %1.4f \n",gradrms);
    if (gradrms<OPTTHRESH) break;
  }

  printf("\n final MM grad RMS: %1.4f \n",gradrms);

//  printf("\n updating IC's \n");
  update_ic();
  for (int i=0;i<nbonds;i++)
  {
    int a1,a2;
    a1 = bonds[i][0];
    a2 = bonds[i][1];
    if (bond_stretch(a1,a2)/ffbondd(a1,a2)>1.5)
      printf(" warning: bond %i %i far from eq. dist %1.2f, actual: %1.2f \n",a1,a2,ffbondd(a1,a2),bond_stretch(a1,a2));
  }

//  printf(" new XYZ \n");
//  print_xyz();
//  print_ic();
//  printf("\n");

  delete [] dxm1;

  if (gradrms>1 || gradrms!=gradrms)
  {
    printf(" opt fails \n");
    return 0;
  }
  else 
    return 1;
}

//Shadow mixes in old IC into gradient
int ICoord::opt(string xyzfile_string, ICoord shadow){

//  printf("  \n"); 
  dxm1 = new double[3*natoms];

  mm_grad(shadow);
  //print_grad();

  ofstream xyzfile;
//  string xyzfile_string = "xyzfile.txt";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

  for (int n=0;n<OPTSTEPS;n++)
  {
//    printf(" Opt step: %i ",n+1);

    xyzfile << " " << natoms << endl;
    xyzfile << " " << gradrms << endl;
    for (int i=0;i<natoms;i++) 
    {
      xyzfile << "  " << anames[i];
      xyzfile << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2];
      xyzfile << endl;
    }

    if (n%100==0) 
      pgradrms = 10000; // reset CG optimizer

    if (n==SHADOWSTEPS) 
      pgradrms = 10000; // reset CG optimizer
    if (n==0) update_xyz_sd();
    else update_xyz_cg();

    if (n<SHADOWSTEPS)
    {
      for (int i=0;i<3*natoms;i++)
        shadow.coords[i] = coords[i];
      mm_grad(shadow);
    }
    else
      mm_grad();
//    print_grad();
//    printf(" gradient RMS: %1.4f \n",gradrms);
    if (gradrms<OPTTHRESH) break;
  }

  printf(" final MM grad RMS: %1.4f \n",gradrms);

//  printf("\n updating IC's \n");
  update_ic();
  for (int i=0;i<nbonds;i++)
  {
    int a1,a2;
    a1 = bonds[i][0];
    a2 = bonds[i][1];
    if (bond_stretch(a1,a2)/ffbondd(a1,a2)>1.5)
      printf(" warning: bond %i %i far from eq. dist %1.2f, actual: %1.2f \n",a1,a2,ffbondd(a1,a2),bond_stretch(a1,a2));
  }

//  printf(" new XYZ \n");
//  print_xyz();
//  print_ic();
//  printf("\n");

  delete [] dxm1;

  if (gradrms>1 || gradrms!=gradrms)
  {
    printf(" opt fails \n");
    return 0;
  }
  else 
    return 1;
}

// steepest descent
void ICoord::update_xyz_sd(){
  
  for (int i=0;i<3*natoms;i++)
    if(abs(grad[i])>MAXAD)
      grad[i]=sign(grad[i])*0.1;
  double SCALE = SCALESD;
  for (int i=0;i<3*natoms;i++)
  {
    dxm1[i] = grad[i]/SCALE;
    coords[i] += dxm1[i];
  }

  return;
} 

void ICoord::update_xyz_cg(){
  
  for (int i=0;i<3*natoms;i++)
    if(abs(grad[i])>MAXAD)
      grad[i]=sign(grad[i])*MAXAD;
  double SCALE = SCALE1; // 1.5 was fine
  double SCALE2 = gradrms*gradrms/pgradrms/pgradrms / SCALEA;
  if (SCALE2 > 1/SCALEA) SCALE2 = 1/SCALEA;
//  printf(" GRADRMS: %1.3f SCALE2: %1.3f \n",gradrms,SCALE2);
  for (int i=0;i<3*natoms;i++)
    coords[i] += grad[i]/SCALE+dxm1[i]*SCALE2;

  return;
} 

