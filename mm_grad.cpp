// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "icoord.h"
#include "utils.h"
using namespace std;

#define SHADOWFACTOR 1.5

// this function will set up MM parameter arrays
void ICoord::mm_init(){

//  printf(" in mm_init, setting default values \n");

  for (int i=0;i<natoms;i++) 
  {
    if (anumbers[i] == 1)
    {
//VDW 1.2
//CHARM: 1.34 A, 0.03 kcal/mol
      ffR[i] = 1.2;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 2)
    {
//VDW 1.4
      ffR[i] = 1.4;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 3)
    {
//VDW 1.4
      ffR[i] = 2.5;
      ffeps[i] = 0.003;
    }
    else if (anumbers[i] == 5)
    {
//VDW N/A, using same as carbon
      ffR[i] = 1.7;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 6)
    {
//VDW 1.7
//CHARM 2.1, 0.07
      ffR[i] = 1.7;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 7)
    {
//VDW 1.55
//CHARM 1.85, 0.2
      ffR[i] = 1.55;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 8)
    {
//VDW 1.52
//CHARMM 1.7, 0.12
      ffR[i] = 1.52;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 9)
    {
//VDW 147
//CHARMM 1.62 A, 0.1 kcal/mol
      ffR[i] = 1.6;
      ffeps[i] = 0.0005;
    }
    else if (anumbers[i] == 11) //Na
    {
      ffR[i] = 3.0;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 12) //Mg
    {
      ffR[i] = 3.0;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 13)
    {
      ffR[i] = 2.6;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 14)
    {
      ffR[i] = 2.1;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 15)
    {
      ffR[i] = 1.8;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 16)
    {
      ffR[i] = 1.8;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 17)
    {
      ffR[i] = 1.8;
      ffeps[i] = 0.002;
    }
    else if (anumbers[i] == 26)
    {
      ffR[i] = 2.5;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 27)
    {
      ffR[i] = 2.5;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 28)
    {
      ffR[i] = 2.6;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 29)
    {
      ffR[i] = 2.7;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 35)
    {
      ffR[i] = 3.7;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 40)
    {
      ffR[i] = 3.5;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 46)
    {
      ffR[i] = 3.0;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 47)
    {
      ffR[i] = 2.9;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 48)
    {
      ffR[i] = 3.0;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 49)
    {
      ffR[i] = 3.1;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 73)
    {
      ffR[i] = 3.3;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 77)
    {
      ffR[i] = 3.5;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 78)
    {
      ffR[i] = 3.5;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 79)
    {
      ffR[i] = 3.6;
      ffeps[i] = 0.005;
    }
    else if (anumbers[i] == 0)
    {
      ffR[i] = 1.5;
      ffeps[i] = 0.00;
    }
    else
      printf(" mm_grad not implemented for element %2i \n",anumbers[i]);
  }

  return;
}

double ICoord::mm_energy()
{
  double E = 0;
  E += vdw_energy_all();
  //others not yet implemented

  return E;
}

int ICoord::mm_grad(){

  for (int i=0;i<3*natoms;i++)
    grad[i] = 0;

//  bond_grad_all();
//  angle_grad_all();
  vdw_grad_all();  

  //torsion_grad not yet complete
  //torsion_grad_all();

//  imptor_grad_all();

  pgradrms = gradrms;
  gradrms = 0;
  for (int i=0;i<3*natoms;i++)
    gradrms += grad[i]*grad[i];
  gradrms = sqrt(gradrms);

  return 0;
}


int ICoord::mm_grad(ICoord shadow){

  for (int i=0;i<3*natoms;i++)
    grad[i] = 0;
  for (int i=0;i<3*natoms;i++)
    shadow.grad[i] = 0;

  bond_grad_all();
  angle_grad_all();
  vdw_grad_all();  

  //torsion_grad not yet complete
  //torsion_grad_all();

  imptor_grad_all();


  shadow.mm_grad();
  for (int i=0;i<3*natoms;i++)
    grad[i] += SHADOWFACTOR*shadow.grad[i];

  pgradrms = gradrms;
  gradrms = 0;
  for (int i=0;i<3*natoms;i++)
    gradrms += grad[i]*grad[i];
  gradrms = sqrt(gradrms);

  return 0;
}


void ICoord::vdw_grad_all(){

  for (int i=0;i<n_nonbond;i++)
    vdw_grad_1(nonbond[i][0],nonbond[i][1]);

  return;
}

void ICoord::bond_grad_all(){

  for (int i=0;i<nbonds;i++)
    bond_grad_1(bonds[i][0],bonds[i][1]);

  return;
}

void ICoord::angle_grad_all(){

  for (int i=0;i<nangles;i++)
    angle_grad_1(angles[i][0],angles[i][1],angles[i][2]);

  return;
}

void ICoord::torsion_grad_all(){

  for (int i=0;i<ntor;i++)
    torsion_grad_1(torsions[i][0],torsions[i][1],torsions[i][2],torsions[i][3]);

  return;
}

void ICoord::imptor_grad_all(){

  for (int i=0;i<nimptor;i++)
    imptor_grad_1(imptor[i][0],imptor[i][1],imptor[i][2],imptor[i][3]);

  return;
}

double ICoord::vdw_energy_all()
{
  double E = 0;
  for (int i=0;i<n_nonbond;i++)
    E += vdw_energy_1(nonbond[i][0],nonbond[i][1]);

  return E;
}

double ICoord::vdw_energy_1(int i, int j)
{
  double R = ffR[i] + ffR[j];
  double eps = sqrt( ffeps[i] * ffeps[j] );

  double r = distance(i,j);

//  printf(" R: %1.4f r: %1.4f \n",R,r);
  double Rr = R / r;
  double Rr6 = Rr*Rr*Rr;
  Rr6 = Rr6*Rr6;

  double E = eps * ( Rr6*Rr6 - 2*Rr6 );

  return E;
}

void ICoord::vdw_grad_1(int i, int j){

  double R = ffR[i] + ffR[j];
  double eps = sqrt( ffeps[i] * ffeps[j] );

  double* dx = new double[3];
  dx[0] = coords[3*i+0]-coords[3*j+0];
  dx[1] = coords[3*i+1]-coords[3*j+1];
  dx[2] = coords[3*i+2]-coords[3*j+2];

  double r = distance(i,j);

//  printf(" R: %1.4f r: %1.4f \n",R,r);
  double Rr = R / r;
  double Rr2 = Rr*Rr;
  double Rr6 = Rr2*Rr2*Rr2;
  
  double t = - eps * Rr6 * Rr2 * (-12 * Rr6 + 12) / R / R;
//  printf(" t on %i %i is: %1.4f \n",i,j,t);
 // t *= -1; //CPMZ TMP
  if (t>0.5)
  {
   // printf(" t: %4.1f",t);
   t = 0.5;
  }
  grad[3*i+0] += t*dx[0];
  grad[3*i+1] += t*dx[1];
  grad[3*i+2] += t*dx[2];
  grad[3*j+0] -= t*dx[0];
  grad[3*j+1] -= t*dx[1];
  grad[3*j+2] -= t*dx[2];
 
  delete [] dx;

  return;
}

double ICoord::bond_stretch(int i, int j) {
 
  return distance(i,j) - ffbondd(i,j);
}

void ICoord::bond_grad_1(int i, int j){

  double* dx = new double[3];
  dx[0] = coords[3*i+0]-coords[3*j+0];
  dx[1] = coords[3*i+1]-coords[3*j+1];
  dx[2] = coords[3*i+2]-coords[3*j+2];

//  double d = ( distance(i,j) - ffbondd(i,j) );
  double d = bond_stretch(i,j);
  
  double t = - d * ffbonde(i,j);
  //printf(" magnitude of gradient on %i %i: %1.4f, current stretch: %1.4f \n",i,j,t,d);
  //printf(" equilibrium distance: %1.4f energy: %1.4f \n",ffbondd(i,j),ffbonde(i,j));
  grad[3*i+0] += t*dx[0];
  grad[3*i+1] += t*dx[1];
  grad[3*i+2] += t*dx[2];
  grad[3*j+0] -= t*dx[0];
  grad[3*j+1] -= t*dx[1];
  grad[3*j+2] -= t*dx[2];
 
  delete [] dx;

  return;
}

void ICoord::angle_grad_1(int i, int j, int k){

  double angle = angle_val(i,j,k) *3.14159/180; // in radians

  double* dx1 = new double[3];
  double* dx2 = new double[3];
  dx1[0] = coords[3*i+0]-coords[3*j+0];
  dx1[1] = coords[3*i+1]-coords[3*j+1];
  dx1[2] = coords[3*i+2]-coords[3*j+2];
  dx2[0] = coords[3*k+0]-coords[3*j+0];
  dx2[1] = coords[3*k+1]-coords[3*j+1];
  dx2[2] = coords[3*k+2]-coords[3*j+2];

  double r1 = distance(i,j); 
  double r2 = distance(j,k);
//  printf(" r1: %1.4f r2: %1.4f\n",r1,r2);
  double* g1 = new double[3];
  double* g2 = new double[3];
  double g1a = -r2/r1*cos(angle);
  double g2a = -r1/r2*cos(angle);
//  printf(" g1a: %1.4f g2a: %1.4f\n",g1a,g2a);
  g1[0] = g1a*dx1[0]+dx2[0];
  g1[1] = g1a*dx1[1]+dx2[1];
  g1[2] = g1a*dx1[2]+dx2[2];
  g2[0] = g2a*dx2[0]+dx1[0];
  g2[1] = g2a*dx2[1]+dx1[1];
  g2[2] = g2a*dx2[2]+dx1[2];

  double d = angle - ffangled(i,j)*3.14/180; //in radians
  
  double SCALE=1;
  double t = 2 * d * ffanglee(i,j) / sin(angle) / r1 / r2 / SCALE;

//  printf(" grad mag on %i %i %i: %1.4f, current stretch: %1.4f angle: %1.4f \n",i,j,k,t,d,angle*180/3.14);
//  printf(" equilibrium distance: %1.4f energy: %1.4f \n",ffangled(i,j),ffanglee(i,j));
  grad[3*i+0] += t*g1[0];
  grad[3*i+1] += t*g1[1];
  grad[3*i+2] += t*g1[2];
  grad[3*j+0] -= t*g1[0];
  grad[3*j+1] -= t*g1[1];
  grad[3*j+2] -= t*g1[2];
  grad[3*j+0] -= t*g2[0];
  grad[3*j+1] -= t*g2[1];
  grad[3*j+2] -= t*g2[2];
  grad[3*k+0] += t*g2[0];
  grad[3*k+1] += t*g2[1];
  grad[3*k+2] += t*g2[2];
 
  delete [] dx1;
  delete [] dx2;
  delete [] g1;
  delete [] g2;

  return;
}

void ICoord::torsion_grad_1(int i, int j, int k, int l){

  double angle = abs(torsion_val(i,j,k,l)*3.14159/180);

//  i
//  jk
//  l
  double* rij = new double[3];
  double* rjk = new double[3];
  double* rkj = new double[3]; 
  double* rlk = new double[3];

  rij[0] = coords[3*i+0] - coords[3*j+0];
  rij[1] = coords[3*i+1] - coords[3*j+1];
  rij[2] = coords[3*i+2] - coords[3*j+2];
  rkj[0] = coords[3*k+0] - coords[3*j+0];
  rkj[1] = coords[3*k+1] - coords[3*j+1];
  rkj[2] = coords[3*k+2] - coords[3*j+2];
  rlk[0] = coords[3*l+0] - coords[3*k+0];
  rlk[1] = coords[3*l+1] - coords[3*k+1];
  rlk[2] = coords[3*l+2] - coords[3*k+2];
  rjk[0] = -rkj[0];
  rjk[1] = -rkj[1];
  rjk[2] = -rkj[2];
  double R2 = rkj[0]*rkj[0]+rkj[1]*rkj[1]+rkj[2]*rkj[2];
  double R = sqrt(R2);
 
  double* m = new double[3];
  double* n = new double[3];

  cross(m,rij,rkj);
  cross(n,rjk,rlk);

  double mm = m[0]*m[0]+m[1]*m[1]+m[2]*m[2];
  double nn = n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
  double dot = rij[0]*n[0]+rij[1]*n[1]+rij[2]*n[2];
  if (dot < 0) angle *= -1;

  double d1 = rij[0]*rkj[0] + rij[1]*rkj[1] + rij[2]*rkj[2];
  double d2 = rlk[0]*rjk[0] + rlk[1]*rjk[1] + rlk[2]*rjk[2];

// this term must be updated:
//fftorm(i,j,k,l)
//  printf(" sin(0): %1.4f \n",sin(angle - fftord(i,j,k,l)*3.14/180));
  double dEdphi = 2 * fftore(i,j,k,l) * sin(angle - fftord(i,j,k,l)*3.14159/180);

  printf(" grad mag on torsion %i %i %i %i: %1.4f, angle: %1.4f \n",i,j,k,l,dEdphi,angle*180/3.14);
  printf(" equilibrium distance: %1.4f energy: %1.4f \n",fftord(i,j,k,l),fftore(i,j,k,l));

  double* S = new double[3];
  double* Fi = new double[3];
  double* Fl = new double[3];
  grad[3*i+0] += Fi[0] = -dEdphi*R*m[0]/mm;
  grad[3*i+1] += Fi[1] = -dEdphi*R*m[1]/mm;
  grad[3*i+2] += Fi[2] = -dEdphi*R*m[2]/mm;
  grad[3*l+0] += Fl[0] = dEdphi*R*n[0]/nn;
  grad[3*l+1] += Fl[1] = dEdphi*R*n[1]/nn;
  grad[3*l+2] += Fl[2] = dEdphi*R*n[2]/nn;
 
  S[0] = (Fi[0]*d1-Fl[0]*d2)/R2;
  S[1] = (Fi[1]*d1-Fl[1]*d2)/R2;
  S[2] = (Fi[2]*d1-Fl[2]*d2)/R2;

  grad[3*j+0] += -Fi[0] + S[0];
  grad[3*j+1] += -Fi[1] + S[1];
  grad[3*j+2] += -Fi[2] + S[2];
  grad[3*k+0] -=  Fl[0] + S[0];
  grad[3*k+1] -=  Fl[1] + S[1];
  grad[3*k+2] -=  Fl[2] + S[2];

  delete [] rij;
  delete [] rjk;
  delete [] rkj;
  delete [] rlk;
  delete [] m;
  delete [] n;
  delete [] S;
  delete [] Fi;
  delete [] Fl;

  return;
}

void ICoord::imptor_grad_1(int i, int j, int k, int l){

//w is old ijkl order
//w  double angle = abs(torsion_val(i,j,k,l)*3.14159/180);
  double angle = abs(torsion_val(j,k,i,l)*3.14159/180);

//  k is central
//  ij
//  k
//  l
  double* rij = new double[3];
  double* rjk = new double[3];
  double* rkj = new double[3]; 
  double* rlk = new double[3];

  rij[0] = coords[3*i+0] - coords[3*j+0];
  rij[1] = coords[3*i+1] - coords[3*j+1];
  rij[2] = coords[3*i+2] - coords[3*j+2];
  rkj[0] = coords[3*k+0] - coords[3*j+0];
  rkj[1] = coords[3*k+1] - coords[3*j+1];
  rkj[2] = coords[3*k+2] - coords[3*j+2];
  rlk[0] = coords[3*l+0] - coords[3*k+0];
  rlk[1] = coords[3*l+1] - coords[3*k+1];
  rlk[2] = coords[3*l+2] - coords[3*k+2];
  rjk[0] = -rkj[0];
  rjk[1] = -rkj[1];
  rjk[2] = -rkj[2];
  double R2 = rkj[0]*rkj[0]+rkj[1]*rkj[1]+rkj[2]*rkj[2];
  double R = sqrt(R2);
 
  double* m = new double[3];
  double* n = new double[3];

  cross(m,rij,rkj);
  cross(n,rjk,rlk);

  double mm = m[0]*m[0]+m[1]*m[1]+m[2]*m[2];
  double nn = n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
  double dot = rij[0]*n[0]+rij[1]*n[1]+rij[2]*n[2];
  if (dot < 0) angle *= -1;

  double d1 = rij[0]*rkj[0] + rij[1]*rkj[1] + rij[2]*rkj[2];
  double d2 = rlk[0]*rjk[0] + rlk[1]*rjk[1] + rlk[2]*rjk[2];
//w  double dEdphi = 2 * ffimptore(i,j,k,l) * angle;
//  double dEdphi = 2 * ffimptore(i,j,k,l) * angle;
//old flat only  double dEdphi = -2 * ffimptore(i,j,k,l) * sin(3.14159 - angle);
  double dAngle;
  double dEdphi; 
  if (angle>=0)
  {
    dAngle = ffimptord(i,j,k,l) - angle;
    dEdphi = -2 * ffimptore(i,j,k,l) * sin(dAngle);
  //  printf(" dAngle1: %1.1f ffimptord: %1.1f angle: %1.1f \n",dAngle,ffimptord(i,j,k,l),angle);
  }
  else
  {
    dAngle = (ffimptord(i,j,k,l) + angle);
    dEdphi = -2 * ffimptore(i,j,k,l) * sin(-dAngle);
  //  printf(" dAngle2: %1.1f ffimptord: %1.1f angle: %1.1f \n",dAngle,-ffimptord(i,j,k,l),angle);
  }

//  printf(" grad mag on imptor %i %i %i %i: %1.4f, angle: %1.4f \n",i,j,k,l,dEdphi,angle*180/3.14);
//  printf(" equilibrium distance: %1.4f energy: %1.4f \n",180.0,ffimptore(i,j,k,l));

  double* S = new double[3];
  double* Fi = new double[3];
  double* Fl = new double[3];
  grad[3*i+0] += Fi[0] = -dEdphi*R*m[0]/mm;
  grad[3*i+1] += Fi[1] = -dEdphi*R*m[1]/mm;
  grad[3*i+2] += Fi[2] = -dEdphi*R*m[2]/mm;
  grad[3*l+0] += Fl[0] = dEdphi*R*n[0]/nn;
  grad[3*l+1] += Fl[1] = dEdphi*R*n[1]/nn;
  grad[3*l+2] += Fl[2] = dEdphi*R*n[2]/nn;
 
  S[0] = (Fi[0]*d1-Fl[0]*d2)/R2;
  S[1] = (Fi[1]*d1-Fl[1]*d2)/R2;
  S[2] = (Fi[2]*d1-Fl[2]*d2)/R2;

  grad[3*j+0] += -Fi[0] + S[0];
  grad[3*j+1] += -Fi[1] + S[1];
  grad[3*j+2] += -Fi[2] + S[2];
  grad[3*k+0] -=  Fl[0] + S[0];
  grad[3*k+1] -=  Fl[1] + S[1];
  grad[3*k+2] -=  Fl[2] + S[2];

  delete [] rij;
  delete [] rjk;
  delete [] rkj;
  delete [] rlk;
  delete [] m;
  delete [] n;
  delete [] S;
  delete [] Fi;
  delete [] Fl;

  return;
}

double ICoord::ffbondd(int i, int j){

  if (anumbers[i]==1 && anumbers[j]==6)
    return 1.09;
  else if (anumbers[i]==6 && anumbers[j]==1)
    return 1.09;
  else if (anumbers[i]==6 && anumbers[j]==6)
    return 1.5;
  else if (anumbers[i]==1 && anumbers[j]==1)
    return 0.8;
  else if (anumbers[i]==1 && anumbers[j]==7)
    return 1.1;
  else if (anumbers[i]==7 && anumbers[j]==1)
    return 1.1;
  else if (anumbers[i]==6 && anumbers[j]==7)
    return 1.45;
  else if (anumbers[i]==7 && anumbers[j]==6)
    return 1.45;
  else if (anumbers[i]==7 && anumbers[j]==7)
    return 1.4;
  else if (anumbers[i]==1 && anumbers[j]==8)
    return 1.0;
  else if (anumbers[i]==8 && anumbers[j]==1)
    return 1.0;
  else if (anumbers[i]==6 && anumbers[j]==8)
    return 1.4;
  else if (anumbers[i]==8 && anumbers[j]==6)
    return 1.4;
  else if (anumbers[i]==1 && anumbers[j]==5)
    return 1.2;
  else if (anumbers[i]==5 && anumbers[j]==1)
    return 1.2;
  else if (anumbers[i]==5 && anumbers[j]==5)
    return 1.5;
  else if (anumbers[i]==6 && anumbers[j]==9)
    return 1.33;
  else if (anumbers[i]==9 && anumbers[j]==6)
    return 1.33;
  else if (anumbers[i]==1 && anumbers[j]==13)
    return 1.6;
  else if (anumbers[i]==13 && anumbers[j]==1)
    return 1.6;
  else if (anumbers[i]==1 && anumbers[j]==15)
    return 1.44;
  else if (anumbers[i]==15 && anumbers[j]==1)
    return 1.44;
  else if (anumbers[i]==1 && anumbers[j]==16)
    return 1.35;
  else if (anumbers[i]==16 && anumbers[j]==1)
    return 1.35;
  else if (anumbers[i]==1 && anumbers[j]==17)
    return 1.33;
  else if (anumbers[i]==17 && anumbers[j]==1)
    return 1.33;
  else if (anumbers[i]==1 && anumbers[j]==73)
    return 1.75;
  else if (anumbers[i]==73 && anumbers[j]==1)
    return 1.75;
  else if (anumbers[i]==7 && anumbers[j]==73)
    return 2.15;
  else if (anumbers[i]==73 && anumbers[j]==7)
    return 2.15;
  else if (anumbers[i]==8 && anumbers[j]==73)
    return 1.9;
  else if (anumbers[i]==73 && anumbers[j]==8)
    return 1.9;
  else
  {
//    printf(" no bond params found for %i %i \n",i,j);
    return (getR(i)+getR(j))/2.2; //was 2
  }
}

double ICoord::ffbonde(int i, int j){
//currently in HT
  return 0.3;
}

double ICoord::ffangled(int i, int j){

  if (anumbers[i]==1 && anumbers[j]==6)
    return 109;
  else if (anumbers[i]==6 && anumbers[j]==1)
    return 109;
  else if (anumbers[i]==6 && anumbers[j]==6)
    return 109;
  else if (anumbers[i]==1 && anumbers[j]==1)
    return 109;
  else
    return 109;
}

double ICoord::ffanglee(int i, int j){

  return 0.05;
}

double ICoord::fftord(int i, int j, int k, int l){

   return 120;
}

double ICoord::fftore(int i, int j, int k, int l){

   return 0.05;
}

double ICoord::fftorm(int i, int j, int k, int l){

   return 3;
}


double ICoord::ffimptore(int i, int j, int k, int l){

//  printf(" in imptore, ijkl: %i %i %i %i \n",i,j,k,l);
//  printf(" k anumber: %i \n",anumbers[k]);
#if 1
  if (anumbers[k]==7)
    return 0.15;
  else
#endif
    return 0.05;
}

double ICoord::ffimptord(int i, int j, int k, int l){

//  printf(" k anumber: %i \n",anumbers[k]);
  if (anumbers[k]==7)
    return 2.1;
  else
    return 3.14159;
}

void ICoord::print_grad(){
  
  printf(" Gradient:\n");
  for (int i=0;i<natoms;i++)
  {
    printf(" %1.4f %1.4f %1.4f \n",grad[3*i+0],grad[3*i+1],grad[3*i+2]);
  }

  return;
}
