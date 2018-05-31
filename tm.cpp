// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "zstruct.h"
#include "constants.h"

#define OVERRIDE_TYPE 0
#define ANGLE_LINEAR 150.0

//1. test 3 code (2 code probably okay)
//2. test 4 code 
//3. write 5 code --> type 0 (need to test alignment)
//4. in add_align_v TM vs not update (done)
//5. double check eta2 breaking tests... ID correct opposites etc


//Note: could call "run" for multiple TM's

void ZStruct::find_linear_pair(int& a1o, int& a2o, int& a1f, int& a2f, int a1, int a2, int a3, int a4, int wtm, ICoord ic1)
{
  double* ang = new double[6];
  ang[0] = ic1.angle_val(a1,wtm,a2);
  ang[1] = ic1.angle_val(a1,wtm,a3);
  ang[2] = ic1.angle_val(a1,wtm,a4);
  ang[3] = ic1.angle_val(a2,wtm,a3);
  ang[4] = ic1.angle_val(a2,wtm,a4);
  ang[5] = ic1.angle_val(a3,wtm,a4);

  int nmax = -1;
  double max = -1.;
  for (int i=0;i<6;i++)
  if (ang[i]>max)
  {
    max = ang[i];
    nmax = i;
  }

  if (nmax==0)
  {
    a1o = a1;
    a2o = a2;
    a1f = a3;
    a2f = a4;
  }
  else if (nmax==1)
  {
    a1o = a1;
    a2o = a3;
    a1f = a2;
    a2f = a4;
  }
  else if (nmax==2)
  {
    a1o = a1;
    a2o = a4;
    a1f = a2;
    a2f = a3;
  }
  else if (nmax==3)
  {
    a1o = a2;
    a2o = a3;
    a1f = a1;
    a2f = a4;
  }
  else if (nmax==4)
  {
    a1o = a2;
    a2o = a4;
    a1f = a1;
    a2f = a3;
  }
  else if (nmax==5)
  {
    a1o = a3;
    a2o = a4;
    a1f = a1;
    a2f = a2;
  }

  delete [] ang;

  return;
}

int ZStruct::run_4c_om(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv)
{
  int nfound = 0;

  printf("  in run_4c_om for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nadd!=2) return 0;
  if (nbrks!=0) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  int** btm = ic1.bondstm;

  int a1 = btm[0][1];
  int a2 = btm[1][1];
  int a3 = btm[2][1];
  int a4 = btm[3][1];
  int a5 = adds[0];
  if (adds[0]==wtm)
    a5 = adds[1];
  int a6 = adds[2];
  if (adds[2]==wtm)
    a6 = adds[3];

  if (nadd==2) //0 or 1 break
  {
    //find linear pair
    int a1o,a2o; //linear pair
    int a1f,a2f; //other two atoms
    find_linear_pair(a1o,a2o,a1f,a2f,a1,a2,a3,a4,wtm,ic1);

    ntor[0] = 0;
    nangles[0] = 2;
    angles[0][0] = a1f;
    angles[0][1] = wtm;
    angles[0][2] = a6;
    angles[1][0] = a2f;
    angles[1][1] = wtm;
    angles[1][2] = a5;

    anglev[0] = anglev[1] = 180.;

    ntor[1] = 0;
    nangles[1] = 2;
    angles[2][0] = a1f;
    angles[2][1] = wtm;
    angles[2][2] = a6;
    angles[3][0] = a2f;
    angles[3][1] = wtm;
    angles[3][2] = a5;

    anglev[2] = anglev[3] = 180.;

    printf("   nadd: 2. angles: %i %i %i and %i %i %i \n",a1f+1,wtm+1,a5+1,a2f+1,wtm+1,a6+1);

    nfound = 2;
    napp[0] = napp[1] = 1; 

    get_vec(wtm,a2f,ic1.coords,&aprv[3]);
    get_vec(wtm,a1f,ic1.coords,&aprv[6]);
  }
  else
    nfound = 0;


  return nfound;
}



int ZStruct::tm_angle_drive(int nr, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* active1)
{
  printf("\n in tm_angle_drive for nadd: %i nbrks: %i  \n",nadd,nbrks);
  int nfound = 0;


  for (int i=0;i<nadd;i++)
    printf("  adding: %2i %2i \n",adds[2*i+0]+1,adds[2*i+1]+1);
  for (int i=0;i<nbrks;i++)
    printf("  brking: %2i %2i \n",brks[2*i+0]+1,brks[2*i+1]+1);
//  ic1.print_xyz();
//  ic1.print_bonds();

 //find up to two metal centers
  int tm1 = -1; int tm2 = -1;
  for (int i=0;i<ic1.natoms;i++)
  if (ic1.isTM(i))
  {
    tm1 = i;
    break;
  }
  for (int i=tm1;i<ic1.natoms;i++)
  if (ic1.isTM(i) && tm1!=i)
  {
    tm2 = i;
    break;
  }


  //temporary bimetallic hack, avoids double TM moves
  int tmtmadd = 0;
 //TM-TM bond add
  for (int i=0;i<nadd;i++)
  if (ic1.isTM(adds[2*i+0]) && ic1.isTM(adds[2*i+1]))
    tmtmadd++;
 //2 different TM adds
  for (int i=0;i<nadd;i++)
  for (int j=0;j<i;j++)
  if (ic1.isTM(adds[2*i+0])+ic1.isTM(adds[2*i+1]) &&
      ic1.isTM(adds[2*j+0])+ic1.isTM(adds[2*j+1]))
  {
    int tm1 = adds[2*i+0];
    if (!ic1.isTM(tm1))
      tm1 = adds[2*i+1];
    int tm2 = adds[2*j+0];
    if (!ic1.isTM(tm2))
      tm2 = adds[2*j+1];
    if (tm1!=tm2)
      tmtmadd++;
  }
  if (tmtmadd)
  {
    printf("   TM-TM bond or 2x TM add, skipping this one \n");
    return 0;
  }

  /* printf("  active:");
  for (int i=0;i<ic1.natoms;i++)
    printf(" %2i",active1[i]);
  printf("\n"); */

  int naddtm = 0;
  int nbrktm = 0;
  int* addtm = new int[2*nadd+1];
  int* brktm = new int[2*nbrks+1];

  int wtm = get_wtm_addbrktm(nadd,adds,nbrks,brks,ic1,naddtm,addtm,nbrktm,brktm);
  //printf("   TM(%s) is atom %i \n",ic1.anames[wtm].c_str(),wtm+1);
  int stm;
  int coordntm = -1;
  int geomtype = -1;
  if (wtm==tm1)
  {
    stm = 1;
    coordntm = ic1.nbondstm;
    geomtype = determine_geomtype(wtm,stm,coordntm,ic1);
  }
  else
  {
    stm = 2;
    coordntm = ic1.nbondstm2;
    geomtype = determine_geomtype(wtm,stm,coordntm,ic1); 
  }
  printf("  geomtype is: %2i for TM %i with coordntm: %i \n",geomtype,wtm+1,coordntm);


  int* napp = new int[10];
  double* aprv = new double[50]; 
  int nroute = run_tm(coordntm,geomtype,wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv,active1);


#if 1
  printf("   found (%i) angle drivers:",nroute);
  for (int nf=0;nf<nroute;nf++)
    printf(" %i",nangles[nf]);
  printf("\n");
  int angc = 0; //count total angles
  for (int nf=0;nf<nroute;nf++)
  {
    for (int i=0;i<nangles[nf];i++)
    {
      printf("    %i:  %2i %2i %2i  %4.1f \n",nf+1,angles[angc][0]+1,angles[angc][1]+1,angles[angc][2]+1,anglev[angc]);
      angc++;
    }
  }
#endif
#if 0
  printf("   found torsion drivers:");
  for (int nf=0;nf<nroute;nf++)
    printf(" %i",ntor[nf]);
  printf("\n");
  int torc = 0; //count total torsions
  for (int nf=0;nf<nroute;nf++)
  for (int i=0;i<ntor[nf];i++)
    printf("    %i:  %2i %2i %2i %2i %4.1f \n",nf+1,tor[torc+i][0]+1,tor[torc+i][1]+1,tor[torc+i][2]+1,tor[torc+i][3],torv[torc++ +i]);
#endif


 //align and create input files
  angc = 0; //torc = 0;
  int nappd = 0;
  for (int i=0;i<nroute;i++)
  {
    if (nr>1)
    for (int j=0;j<napp[i];j++)
    {
      printf("   writing input #%i                       (niso: %4i) \n",j+1,niso);
      int index = 3*(nappd+1); nappd++;
      int tm_align_used = align1.add_align_v(nadd,adds,wtm,&aprv[index]);
      write_initial_xyz(niso,ic1.natoms,ic1.anames,align1.xyza,ic1.q1);
      write_ISOMER(niso++,nadd,adds,nbrks,brks,nangles[i],&angles[angc],&anglev[angc],0,NULL,NULL);
      //exit(1);
      nfound++;
      if (!tm_align_used) break;
    }
    else
    {
      printf("   writing input #%i                       (niso: %4i) \n",i+1,niso);
      write_initial_xyz(niso,ic1.natoms,ic1.anames,ic1.coords,ic1.q1);
      write_ISOMER(niso++,nadd,adds,nbrks,brks,nangles[i],&angles[angc],&anglev[angc],0,NULL,NULL);
      nfound++;
    }
    angc += nangles[i];
    //torc += ntor[i];
  }


  delete [] addtm;
  delete [] brktm;
  delete [] napp;
  delete [] aprv;

  return nfound;
}

int ZStruct::run_tm(int coordntm, int geomtype, int wtm, int stm, int naddtm, int* addtm, int nbrktm, int* brktm, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1)
{
  int nroute = 0;

#if OVERRIDE_TYPE
  nroute = run_4c_om(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv);
#else
  if (coordntm==2) 
  {
    if (geomtype==1)
      nroute = run_2c_1(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv);
    else
      nroute = run_2c_0(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv);
  }
  else if (coordntm==3) 
  {
    if (geomtype==1)
      nroute = run_3c_1(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv);
    else
      nroute = run_3c_0(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv);
  }
  else if (coordntm==4) 
  {
    if (geomtype==1)
      nroute = run_4c_1(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv,active1);
    else
      nroute = run_4c_0(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv,active1);
  }
  else if (coordntm==5) 
  {
    if (geomtype==1)
      nroute = run_5c_1(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv,active1);
    else
    {
      nroute = run_5c_0(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv,active1);
    }
  }
  else if (coordntm==6) 
  {
    if (geomtype==1)
      nroute = run_6c_1(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv,active1);
    else
    {
      printf("   coordntm 6 geomtype 0 not implemented, using type 1 instead \n");
      nroute = run_6c_1(wtm,stm,naddtm,addtm,nbrktm,brktm,ic1,nangles,angles,anglev,ntor,tor,torv,napp,aprv,active1);
//      exit(1);
    }
  }
  else
    printf("  haven't yet implemented coordntm: %i \n",coordntm);
#endif

  return nroute;
}

int ZStruct::run_6c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1)
{
  int nfound = 0;

  printf("  in run_6c_1 for naddtm: %i nbrktm: %i \n",nadd,nbrks);


  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;
  if (nadd>0 && nbrks==0) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;
  int a1 = btm[0][1];
  int a2 = btm[1][1];
  int a3 = btm[2][1];
  int a4 = btm[3][1];
  int a5 = btm[4][1];
  int a6 = btm[5][1];

  //int p11,p12,p13,p14; //avoided using this
  //int p21,p22,p23,p24;
  //int p31,p32,p33,p34;
  //find_oct_planes(wtm,ic1,p11,p12,p13,p14,p21,p22,p23,p24,p31,p32,p33,p34);

  if (eta2_double_break(wtm,stm,6,nbrks,brks,ic1))
    nbrks--;

  if (nadd==1 && nbrks==1) 
  {
    int a7 = adds[0];
    if (adds[0]==wtm)
      a7 = adds[1];
    int b1 = brks[0];
    if (brks[0]==wtm)
      b1 = brks[1];

    int b1o = find_opp(b1,wtm,stm,ic1);
    int b1s1,b1s2;
    find_2_at_90(b1s1,b1s2,b1,b1o,wtm,stm,ic1);
    printf("   opposite atom: %2i(%s) 90's: %i %i \n",b1o+1,ic1.anames[b1o].c_str(),b1s1+1,b1s2+1);

    ntor[0] = 0;
    nangles[0] = 1;

    angles[0][0] = a7;
    angles[0][1] = wtm;
    angles[0][2] = b1o;

    anglev[0] = 180.;

    nfound = 1;
    napp[0] = napp[1] = 4;
    get_vec_45p(b1,wtm,b1s1,ic1.coords,&aprv[3]);
    get_vec_45m(b1,wtm,b1s1,ic1.coords,&aprv[6]);
    get_vec_45p(b1,wtm,b1s2,ic1.coords,&aprv[9]);
    get_vec_45m(b1,wtm,b1s2,ic1.coords,&aprv[12]);

  } //add1/brk1
  else if (nadd==1 && nbrks==2) 
  {
    int a7 = adds[0];
    if (adds[0]==wtm)
      a7 = adds[1];
    int b1 = brks[0];
    if (brks[0]==wtm)
      b1 = brks[1];
    int b2 = brks[2];
    if (brks[2]==wtm)
      b2 = brks[3];

    int b1o = find_opp(b1,wtm,stm,ic1);
    int b2o = find_opp(b2,wtm,stm,ic1);
    int b1s1,b1s2,b2s1,b2s2;
    find_2_at_90(b1s1,b1s2,b1,b1o,wtm,stm,ic1);
    find_2_at_90(b2s1,b2s2,b2,b2o,wtm,stm,ic1);
    printf("   opposite atom: %2i(%s) 90's: %i %i \n",b1o+1,ic1.anames[b1o].c_str(),b1s1+1,b1s2+1);
    printf("   opposite atom: %2i(%s) 90's: %i %i \n",b2o+1,ic1.anames[b2o].c_str(),b2s1+1,b2s2+1);

    ntor[0] = ntor[1] = 0;
    nangles[0] = nangles[1] = 1;

    angles[0][0] = a7;
    angles[0][1] = wtm;
    angles[0][2] = b1o;

    anglev[0] = 180.;

    angles[1][0] = a7;
    angles[1][1] = wtm;
    angles[1][2] = b2o;

    anglev[1] = 180.;

    nfound = 2;
    napp[0] = napp[1] = 4;
    get_vec_45p(b1,wtm,b1s1,ic1.coords,&aprv[3]);
    get_vec_45m(b1,wtm,b1s1,ic1.coords,&aprv[6]);
    get_vec_45p(b1,wtm,b1s2,ic1.coords,&aprv[9]);
    get_vec_45m(b1,wtm,b1s2,ic1.coords,&aprv[12]);
    get_vec_45p(b2,wtm,b2s1,ic1.coords,&aprv[15]);
    get_vec_45m(b2,wtm,b2s1,ic1.coords,&aprv[18]);
    get_vec_45p(b2,wtm,b2s2,ic1.coords,&aprv[21]);
    get_vec_45m(b2,wtm,b2s2,ic1.coords,&aprv[24]);

  } //add1/brk2
  else
    nfound = 1;


  return nfound;
}

void ZStruct::find_2_at_90(int& a1s1, int& a1s2, int a1, int a1o, int wtm, int stm, ICoord ic1)
{
  int* attached = new int[12];
  int nattached = find_attached(wtm,stm,ic1,attached);
 
  int found = 0;
  for (int i=0;i<nattached;i++)
  if (attached[i]!=a1 && attached[i]!=a1o && !found)
  {
    double ang1 = ic1.angle_val(a1,wtm,attached[i]);
    if (ang1<ANGLE_LINEAR)
    for (int j=0;j<i;j++)
    if (attached[j]!=a1 && attached[j]!=a1o)
    {
      double ang2a = ic1.angle_val(a1,wtm,attached[j]);
      double ang2b = ic1.angle_val(attached[i],wtm,attached[j]);
      if (ang2a<ANGLE_LINEAR && ang2b<ANGLE_LINEAR)
      {
        a1s1 = attached[i];
        a1s2 = attached[j];
        break;
      }
    }
  }

  delete [] attached;  

  return;
}


int ZStruct::find_opp(int a1, int wtm, int stm, ICoord ic1)
{
  int a1o = 0;

  int* attached = new int[12];
  int nattached = find_attached(wtm,stm,ic1,attached);
 
  double maxang = 0.;
  for (int i=0;i<nattached;i++)
  if (attached[i]!=a1)
  {
    double ang1 = ic1.angle_val(a1,wtm,attached[i]);
    if (ang1>maxang)
    {
      a1o = attached[i];
      maxang = ang1;
    }
  }

  delete [] attached;

  return a1o;
}

void ZStruct::find_tb_5c(int* a, int wtm, int stm, ICoord ic1)
{
  //printf("   in find_tb_5c for wtm %i stm %i \n",wtm,stm);

  int* attached = new int[10];
  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;   
  int nattached = 5;
  for (int i=0;i<5;i++)
    attached[i] = btm[i][1];

  int a1 = -1;
  int a2 = -1;
  double maxang = 0.;
  for (int i=0;i<nattached;i++)
  for (int j=0;j<i;j++)
  {
    double val1 = ic1.angle_val(attached[i],wtm,attached[j]);
    if (val1>maxang)
    {
      a1 = attached[i]; 
      a2 = attached[j];
      maxang = val1;
    }
  }
  a[0] = a1;
  a[1] = a2;

  int nf = 2;
  for (int i=0;i<nattached;i++)
  if (a[0]!=attached[i] && a[1]!=attached[i])
    a[nf++] = attached[i];

  printf("   maximum angle is: %6.3f for %i-%i \n",maxang,a[0]+1,a[1]+1);
  if (nf!=5)
  {
    printf(" ERROR: couldn't determine trigonal planar attachment structure \n");
    exit(1);
  }

  delete [] attached;

  return;
}

int ZStruct::run_5c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1)
{
  int nfound = 0;

  printf("  in run_5c_0 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;
  if (nadd==2 && nbrks==0) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;
  int a1 = btm[0][1];
  int a2 = btm[1][1];
  int a3 = btm[2][1];
  int a4 = btm[3][1];
  int a5 = btm[4][1];

  int* a1tb = new int[5];
  find_tb_5c(a1tb,wtm,stm,ic1);
  int a1t  = a1tb[0]; //above
  int a1b  = a1tb[1]; //below
  int a1s1 = a1tb[2];
  int a1s2 = a1tb[3];
  int a1s3 = a1tb[4];
  delete [] a1tb;
  printf("   a1t a1b: %i %i  a1s's: %i %i %i \n",a1t+1,a1b+1,a1s1+1,a1s2+1,a1s3+1);
  //printf("    a's: %i   %i %i %i %i ao's: %i %i %i %i \n",a1t+1,a2s+1,a3s+1,a4s+1,a5s+1,a2o+1,a3o+1,a4o+1,a5o+1);

  if (eta2_double_break(wtm,stm,5,nbrks,brks,ic1))
    nbrks--;

  if (nadd==1 && nbrks!=1) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];

    ntor[0] = ntor[1] = ntor[2] = 0;
    nangles[0] = nangles[1] = nangles[2] = 1;

   //3 destinations/approaches from between trigonal ligands
    angles[0][0] = a1s1;
    angles[0][1] = wtm;
    angles[0][2] = a6;

    anglev[0] = 180.;

    angles[1][0] = a1s2;
    angles[1][1] = wtm;
    angles[1][2] = a6;

    anglev[1] = 180.;

    angles[2][0] = a1s3;
    angles[2][1] = wtm;
    angles[2][2] = a6;

    anglev[2] = 180.;

    nfound = 3;
    napp[0] = napp[1] = napp[2] = 1;
    get_vec(wtm,a1s1,ic1.coords,&aprv[3]);
    get_vec(wtm,a1s2,ic1.coords,&aprv[6]);
    get_vec(wtm,a1s3,ic1.coords,&aprv[9]);
  }
  else if (nadd==1 && nbrks==1) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];

    int a2l = brks[0];
    if (brks[0]==wtm)
      a2l = brks[1];

   //put in trigonal plane where leaving ligand was
    if (!(a2l==a1t || a2l==a1b))
    {
      int a1s1o = a1s1;
      int a1s2o = a1s2;
      if (a2l==a1s1o)
        a1s1o = a1s3;
      if (a2l==a1s2o)
        a1s2o = a1s3;

      ntor[0] = 0;
      nangles[0] = 2;

      angles[0][0] = a1s1o;
      angles[0][1] = wtm;
      angles[0][2] = a6;

      anglev[0] = 120.;

      angles[1][0] = a1s2o;
      angles[1][1] = wtm;
      angles[1][2] = a6;

      anglev[1] = 120.;

      nfound = 1;
      napp[0] = napp[1] = 2;

      get_vec(wtm,a1s1o,ic1.coords,&aprv[3]);
      get_vec(wtm,a1s2o,ic1.coords,&aprv[6]);
    }
   //swap with axial
    else
    {
      int a1o = a1t;
      if (a2l==a1o)
        a1o = a1b;

      ntor[0] = 0;
      nangles[0] = 1;

      angles[0][0] = a1o;
      angles[0][1] = wtm;
      angles[0][2] = a6;

      anglev[0] = 180.;

      nfound = 1;
      napp[0] = 3;

      get_vec_45p_align(a1s1,wtm,a1s2,a2l,ic1.coords,&aprv[3]);
      get_vec_45p_align(a1s2,wtm,a1s1,a2l,ic1.coords,&aprv[6]);
      get_vec_45p_align(a1s3,wtm,a1s2,a2l,ic1.coords,&aprv[9]);
    }
  }
  else if (nadd==2 && nbrks==1) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];
    int a7 = adds[2];
    if (adds[2]==wtm)
      a7 = adds[3];

    int a2l = brks[0];
    if (brks[0]==wtm)
      a2l = brks[1];

   //leaving group is in trigonal plane
    if (!(a2l==a1t || a2l==a1b))
    {
      int a1s1o = a1s1;
      int a1s2o = a1s2;
      if (a2l==a1s1o)
        a1s1o = a1s3;
      if (a2l==a1s2o)
        a1s2o = a1s3;

      ntor[0] = ntor[1] = 0;
      nangles[0] = nangles[1] = 2;

      angles[0][0] = a1s1o;
      angles[0][1] = wtm;
      angles[0][2] = a6;

      anglev[0] = 180.;

      angles[1][0] = a1s2o;
      angles[1][1] = wtm;
      angles[1][2] = a7;

      anglev[1] = 180.;

     //second attach configuration
      angles[2][0] = a1s1o;
      angles[2][1] = wtm;
      angles[2][2] = a7;

      anglev[2] = 180.;

      angles[3][0] = a1s2o;
      angles[3][1] = wtm;
      angles[3][2] = a6;

      anglev[3] = 180.;

      nfound = 2;
      napp[0] = napp[1] = 2;

      get_vec(wtm,a1s1o,ic1.coords,&aprv[3]);
      get_vec(wtm,a1s2o,ic1.coords,&aprv[6]);
      aprv[9] = aprv[3]; aprv[10] = aprv[4]; aprv[11] = aprv[5];
      aprv[12] = aprv[6]; aprv[13] = aprv[7]; aprv[14] = aprv[8];
    }
   //leaving group is axial
    else
    {
      int a1o = a1t; //group opposite leaving group
      if (a2l==a1o)
        a1o = a1b;

      ntor[0] = ntor[1] = 0;
      nangles[0] = nangles[1] = 2;

      angles[0][0] = a1o;
      angles[0][1] = wtm;
      angles[0][2] = a6;

      anglev[0] = 180.;

      angles[1][0] = a1o;
      angles[1][1] = wtm;
      angles[1][2] = a7;

      anglev[1] = 90.;

      angles[2][0] = a1o;
      angles[2][1] = wtm;
      angles[2][2] = a7;

      anglev[2] = 180.;

      angles[3][0] = a1o;
      angles[3][1] = wtm;
      angles[3][2] = a6;

      anglev[3] = 90.;

      nfound = 2;
      napp[0] = napp[1] = 3;

      get_vec_45p_align(a1s1,wtm,a1s2,a2l,ic1.coords,&aprv[3]);
      get_vec_45p_align(a1s2,wtm,a1s1,a2l,ic1.coords,&aprv[6]);
      get_vec_45p_align(a1s3,wtm,a1s2,a2l,ic1.coords,&aprv[9]);
      aprv[12] = aprv[3]; aprv[13] = aprv[4]; aprv[14] = aprv[5];
      aprv[15] = aprv[6]; aprv[16] = aprv[7]; aprv[17] = aprv[8];
      aprv[18] = aprv[9]; aprv[19] = aprv[10]; aprv[20] = aprv[11];
    }
  }
  else if (nadd==1 && nbrks==2) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];

    int a2l = brks[0];
    if (brks[0]==wtm)
      a2l = brks[1];
    int a3l = brks[2];
    if (brks[2]==wtm)
      a3l = brks[3];

    printf(" INCOMPLETE: nadd==1 nbrks==2 \n");
    printf("  shouldn't be here! \n");
    exit(1);
  } //1 add, 2 brk
 //for nadd==0
  else
    nfound = 1;


  return nfound;
} 



int ZStruct::run_5c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1)
{
  int nfound = 0;

  printf("  in run_5c_1 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;
  if (nadd==2 && nbrks==0) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;
  int a1 = btm[0][1];
  int a2 = btm[1][1];
  int a3 = btm[2][1];
  int a4 = btm[3][1];
  int a5 = btm[4][1];

  int a1t = find_on_top_5c(wtm,stm,ic1);
  int a2s = a1; //in plane
  int a3s = a2;
  int a4s = a3;
  int a5s = a4;
  if (a1==a1t) a2s = a5;
  if (a2==a1t) a3s = a5;
  if (a3==a1t) a4s = a5;
  if (a4==a1t) a5s = a5;
  int a2o, a3o, a4o, a5o;
  find_linear_pairs(a2o,a3o,a4o,a5o,a2s,a3s,a4s,a5s,wtm,ic1);
  //printf("    a's: %i   %i %i %i %i ao's: %i %i %i %i \n",a1t+1,a2s+1,a3s+1,a4s+1,a5s+1,a2o+1,a3o+1,a4o+1,a5o+1);

  if (eta2_double_break(wtm,stm,5,nbrks,brks,ic1))
    nbrks--;

  if (nadd>0 && nbrks==0) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];

    ntor[0] = 0;
    nangles[0] = 1;

    angles[0][0] = a1t;
    angles[0][1] = wtm;
    angles[0][2] = a6;

    anglev[0] = 180.;
    if (nadd==2) anglev[0] = 170.;

    nfound = 1;
    napp[0] = 1;
    get_vec(wtm,a1t,ic1.coords,&aprv[3]);

    int nf = 0;
    int nang = 0;

   //4 approaches with shifted equitorial ligands
    if (active1[a2s] && nadd==1)
    {
      int a2c = a3s;
      if (a3s==a2o)
        a2c = a4s;

      ntor[nf+1] = 0;
      nangles[nf+1] = 2;

      angles[nang+1][0] = a2o;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a6;
      angles[nang+2][0] = a2s;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a1t;

      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;

      napp[nf+1] = 1;
      get_vec_45p_align(a2s,wtm,a2c,a1t,ic1.coords,&aprv[(nf+2)*3]);
      //align_vec(a1t,wtm,ic1.coords,&aprv[(nf+2)*3]);
      nf++; 
      nang += 2;
      nfound++;
    }
    if (active1[a3s] && nadd==1)
    {
      int a2c = a2s;
      if (a3s==a3o)
        a2c = a4s;

      ntor[nf+1] = 0;
      nangles[nf+1] = 2;

      angles[nang+1][0] = a3o;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a6;
      angles[nang+2][0] = a3s;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a1t;

      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;

      napp[nf+1] = 1;
      get_vec_45p_align(a3s,wtm,a2c,a1t,ic1.coords,&aprv[(nf+2)*3]);
      //align_vec(a1t,wtm,ic1.coords,&aprv[(nf+2)*3]);
      nf++; 
      nang += 2;
      nfound++;
    }
    if (active1[a4s] && nadd==1)
    {
      int a2c = a2s;
      if (a2s==a4o)
        a2c = a3s;

      ntor[nf+1] = 0;
      nangles[nf+1] = 2;

      angles[nang+1][0] = a4o;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a6;
      angles[nang+2][0] = a4s;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a1t;

      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;

      napp[nf+1] = 1;
      get_vec_45p_align(a4s,wtm,a2c,a1t,ic1.coords,&aprv[(nf+2)*3]);
      //align_vec(a1t,wtm,ic1.coords,&aprv[(nf+2)*3]);
      nf++; 
      nang += 2;
      nfound++;
    }
    if (active1[a5s] && nadd==1)
    {
      int a2c = a2s;
      if (a2s==a5o)
        a2c = a3s;

      ntor[nf+1] = 0;
      nangles[nf+1] = 2;

      angles[nang+1][0] = a5o;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a6;
      angles[nang+2][0] = a5s;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a1t;

      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;

      napp[nf+1] = 1;
      get_vec_45p_align(a5s,wtm,a2c,a1t,ic1.coords,&aprv[(nf+2)*3]);
      //align_vec(a1t,wtm,ic1.coords,&aprv[(nf+2)*3]);
      nf++; 
      nang += 2;
      nfound++;
    }

  }
  else if (nadd==1 && nbrks==1) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];

    int a2l = brks[0];
    if (brks[0]==wtm)
      a2l = brks[1];

    int a2lo = a2o;
    if (a2l==a3s) a2lo = a3o;
    if (a2l==a4s) a2lo = a4o;
    if (a2l==a5s) a2lo = a5o;

    int a2c = a2s; //perp. in plane, remains
    if (a2s==a2l || a2s==a2lo)
    {
      a2c = a3s;
      if (a3s==a2l || a3s==a2lo)
        a2c = a4s;
    }

    ntor[0] = 0;
    nangles[0] = 1;

   //bottom approach
    angles[0][0] = a1t;
    angles[0][1] = wtm;
    angles[0][2] = a6;

    anglev[0] = 180.;
 
    napp[0] = 1;
    get_vec(wtm,a1t,ic1.coords,&aprv[3]);

    nfound = 1;

   //approach from 45 above/below; move 180 to opposite group
    if (a1t!=a2l)
    {
      ntor[1] = 1;
      nangles[1] = 1;

      angles[1][0] = a2lo;
      angles[1][1] = wtm;
      angles[1][2] = a6;

      anglev[1] = 180.;
 
      napp[1] = 1;
      get_vec_45p_align(a2l,wtm,a2c,a1t,ic1.coords,&aprv[6]);
      //get_vec_45m(a2l,wtm,a2c,ic1.coords,&aprv[9]);

      nfound = 2;
    } //if not swapping with top group

  }
  else if (nadd==2 && nbrks==1) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];
    int a7 = adds[2];
    if (adds[2]==wtm)
      a7 = adds[3];

    int a2l = brks[0];
    if (brks[0]==wtm)
      a2l = brks[1];

    int a2lo = a2o;
    if (a2l==a3s) a2lo = a3o;
    if (a2l==a4s) a2lo = a4o;
    if (a2l==a5s) a2lo = a5o;

    int a2c = a2s; //perp. in plane, remains
    if (a2s==a2l || a2s==a2lo)
    {
      a2c = a3s;
      if (a3s==a2l || a3s==a2lo)
        a2c = a4s;
    }

   //simple approach from bottom
    if (a1t==a2l)
    {
      ntor[0] = 0;
      nangles[0] = 0;

      napp[0] = 1;
      get_vec(wtm,a1t,ic1.coords,&aprv[3]);

      nfound = 1;
    }
   //displace approach
    else
    {
      ntor[0] = ntor[1] = 0;
      nangles[0] = nangles[1] = 2;

      angles[0][0] = a1t;
      angles[0][1] = wtm;
      angles[0][2] = a6;
      angles[1][0] = a2lo;
      angles[1][1] = wtm;
      angles[1][2] = a7;

      anglev[0] = 180.;
      anglev[1] = 180.;

      angles[2][0] = a1t;
      angles[2][1] = wtm;
      angles[2][2] = a7;
      angles[3][0] = a2lo;
      angles[3][1] = wtm;
      angles[3][2] = a6;

      anglev[2] = 180.;
      anglev[3] = 180.;

      napp[0] = napp[1] = 1.;
      get_vec_45m_align(a2l,wtm,a2c,a1t,ic1.coords,&aprv[3]);
      //align_vec(a1t,wtm,ic1.coords,&aprv[3]);
      aprv[6] = aprv[3]; aprv[7] = aprv[4]; aprv[8] = aprv[5];

      nfound = 2;
    }
  }
  else if (nadd==1 && nbrks==2) 
  {
    int a6 = adds[0];
    if (adds[0]==wtm)
      a6 = adds[1];

    int a2l = brks[0];
    if (brks[0]==wtm)
      a2l = brks[1];
    int a3l = brks[2];
    if (brks[2]==wtm)
      a3l = brks[3];

   //approach from bottom
    if (a1t==a2l || a1t==a3l)
    {
      int a2lo; //opposite of leaving group
      if (a1t==a3l)
      {
        a2lo = a2o;
        if (a2l==a3s) a2lo = a3o;
        if (a2l==a4s) a2lo = a4o;
        if (a2l==a5s) a2lo = a5o;
      }
      else
      {
        a2lo = a2o;
        if (a3l==a3s) a2lo = a3o;
        if (a3l==a4s) a2lo = a4o;
        if (a3l==a5s) a2lo = a5o;
      }

      ntor[0] = 0;
      nangles[0] = 1;

      angles[0][0] = a2lo;
      angles[0][1] = wtm;
      angles[0][2] = a6;

      anglev[0] = 180.;

      napp[0] = 1;
      get_vec(wtm,a1t,ic1.coords,&aprv[3]);

      nfound = 1;
    } //if a1t is leaving
    else
    {
      int a2lo = a2o; //opposite of leaving group
      if (a2l==a3s) a2lo = a3o;
      if (a2l==a4s) a2lo = a4o;
      if (a2l==a5s) a2lo = a5o;
      int a3lo = a2o;
      if (a3l==a3s) a2lo = a3o;
      if (a3l==a4s) a2lo = a4o;
      if (a3l==a5s) a2lo = a5o;

      ntor[0] = ntor[1] = 0;
      nangles[0] = nangles[1] = 1;

      angles[0][0] = a2lo;
      angles[0][1] = wtm;
      angles[0][2] = a6;

      anglev[0] = 180.;

      angles[1][0] = a3lo;
      angles[1][1] = wtm;
      angles[1][2] = a6;

      anglev[1] = 180.;

      napp[0] = napp[1] = 1;
      get_vec(wtm,a1t,ic1.coords,&aprv[3]);
      aprv[6] = aprv[3]; aprv[7] = aprv[4]; aprv[8] = aprv[5];

      nfound = 2;
    }
  } //1 add, 2 brk
  else
    nfound = 1;


  return nfound;
} 

int ZStruct::run_4c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1)
{
  int nfound = 0;

  printf("  in run_4c_0 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;
  int a1 = btm[0][1];
  int a2 = btm[1][1];
  int a3 = btm[2][1];
  int a4 = btm[3][1];

  if (eta2_double_break(wtm,stm,4,nbrks,brks,ic1))
    nbrks--;

  if (nadd==1) 
  {
    ntor[0] = 0;
    nangles[0] = 0;

    nfound = 1;
    napp[0] = 4;
    get_vec(a1,wtm,ic1.coords,&aprv[3]);
    get_vec(a2,wtm,ic1.coords,&aprv[6]);
    get_vec(a3,wtm,ic1.coords,&aprv[9]);
    get_vec(a4,wtm,ic1.coords,&aprv[12]);
  }
  else if (nadd==2) 
  {
    int a5 = adds[0];
    if (adds[0]==wtm)
      a5 = adds[1];
    int a6 = adds[2];
    if (adds[2]==wtm)
      a6 = adds[3];

    int* a = new int[4];
    a[0] = a1;
    a[1] = a2;
    a[2] = a3;
    a[3] = a4;
#if 0
    if (nbrks==1)
    {
      if (brks[0]==a1 || brks[1]==a1)
        a[0] = -1;
      if (brks[0]==a2 || brks[1]==a2)
        a[1] = -1;
      if (brks[0]==a3 || brks[1]==a3)
        a[2] = -1;
      if (brks[0]==a4 || brks[1]==a4)
        a[3] = -1;
    }
#endif

    for (int j=0;j<4;j++)
    {
      int i1 = 4*j+0; int i2 = 4*j+1;
      int i3 = 4*j+2; int i4 = 4*j+3;

      ntor[2*j+0] = ntor[2*j+1] = 0;
      nangles[2*j+0] = nangles[2*j+1] = 2;

      angles[i1][0] = a5;
      angles[i1][1] = wtm;
      angles[i1][2] = a6;
      angles[i2][0] = a5;
      angles[i2][1] = wtm;
      angles[i2][2] = a[j];

      anglev[i1] = 90.;
      anglev[i2] = 180.;
  
      angles[i3][0] = a5;
      angles[i3][1] = wtm;
      angles[i3][2] = a6;
      angles[i4][0] = a6;
      angles[i4][1] = wtm;
      angles[i4][2] = a[j];

      anglev[i3] = 90.;
      anglev[i4] = 180.;
    }

    delete [] a;

    nfound = 8;
    napp[0] = napp[1] = napp[2] = napp[3] = napp[4] = napp[5] = napp[6] = napp[7] = 1;
    get_vec(a1,wtm,ic1.coords,&aprv[3]);
    get_vec(a1,wtm,ic1.coords,&aprv[6]);
    get_vec(a2,wtm,ic1.coords,&aprv[9]);
    get_vec(a2,wtm,ic1.coords,&aprv[12]);
    get_vec(a3,wtm,ic1.coords,&aprv[15]);
    get_vec(a3,wtm,ic1.coords,&aprv[18]);
    get_vec(a4,wtm,ic1.coords,&aprv[21]);
    get_vec(a4,wtm,ic1.coords,&aprv[24]);
  }
  else
    nfound = 1;


  return nfound;
}

int ZStruct::run_4c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv, int* active1)
{
  int nfound = 0;

  printf("  in run_4c_1 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  int** btm = ic1.bondstm;
  if (stm==2) 
    btm = ic1.bondstm2;
  int a1 = btm[0][1];
  int a2 = btm[1][1];
  int a3 = btm[2][1];
  int a4 = btm[3][1];

  if (eta2_double_break(wtm,stm,4,nbrks,brks,ic1))
    nbrks--;

  if (nadd==1 && nbrks==0) 
  {
    int a5 = adds[0];
    if (adds[0]==wtm)
      a5 = adds[1];
 
    int a1o, a2o, a3o, a4o;
    find_linear_pairs(a1o,a2o,a3o,a4o,a1,a2,a3,a4,wtm,ic1);

    int a1p = a2;
    if (a2==a1o)
      a1p = a3;

    ntor[0] = 0;
    nangles[0] = 1;

   //add on top
    angles[0][0] = a1;
    angles[0][1] = wtm;
    angles[0][2] = a5;
 
    anglev[0] = 90.;

    nfound = 1;
    napp[0] = 2;
    get_oop_vec(a1,wtm,a1p,ic1.coords,&aprv[3]);
    aprv[6] = -aprv[3]; aprv[7] = -aprv[4]; aprv[8]=-aprv[5];

   //add on side with push
    int nang = 1;
    if (active1[a1o] && 0) //CPMZ incomplete
    {
      ntor[nang] = 0;
      nangles[nang] = 2;

      angles[nang][0]   = a1;
      angles[nang][1]   = wtm;
      angles[nang][2]   = a5;
      angles[nang+1][0] = a1;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a1o;

      anglev[nang]   = 180.;
      anglev[nang+1] = 90.;
      nang += 2;
      nfound++;
    }

    napp[1] = napp[2] = napp[3] = napp[4] = 2; 
   //need to code these approaches
  }
  else if (nadd==1 && nbrks==1) 
  {
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1])
      a3 = btm[2][2];
    if (btm[3][2]==brks[0] || btm[3][2]==brks[1])
      a4 = btm[3][2];
    int a5 = adds[0];
    if (adds[0]==wtm)
      a5 = adds[1];
 
    int a1o, a2o, a3o, a4o;
    find_linear_pairs(a1o,a2o,a3o,a4o,a1,a2,a3,a4,wtm,ic1);

    int a0l = a1;
    int a0lo = a1o; //a0lo is atom opposite leaving atom
    if (brks[0]==a2 || brks[1]==a2)
    {
      a0l = a2;
      a0lo = a2o;
    }
    if (brks[0]==a3 || brks[1]==a3)
    {
      a0l = a3;
      a0lo = a3o;
    }
    if (brks[0]==a4 || brks[1]==a4)
    {
      a0l = a4;
      a0lo = a4o;
    }
    //printf("   a0l: %i a0lo: %i \n",a0l+1,a0lo+1);
    //printf("   a1: %i a2: %i a3: %i a4: %i \n",a1+1,a2+1,a3+1,a4+1);

    int a1c = a1; //remaining atom at 90 degrees
    if (a2!=a0l && a2!=a0lo)
      a1c = a2;
    if (a3!=a0l && a3!=a0lo)
      a1c = a3;
    if (a4!=a0l && a4!=a0lo)
      a1c = a4;

    ntor[0] = 0;
    nangles[0] = 1;

   //add opposite linear pair
    angles[0][0] = a0lo;
    angles[0][1] = wtm;
    angles[0][2] = a5;

    anglev[0] = 180.;

    nfound = 1;
    napp[0] = 2;
    //printf("   a0l: %i wtm: %i a1c: %i \n",a0l+1,wtm+1,a1c+1);
    get_vec_45p(a0l,wtm,a1c,ic1.coords,&aprv[3]);
    get_vec_45m(a0l,wtm,a1c,ic1.coords,&aprv[6]);
  }
  else if (nadd==2 && nbrks==0) 
  {
    int a5 = adds[0];
    if (adds[0]==wtm)
      a5 = adds[1];
    int a6 = adds[2];
    if (adds[2]==wtm)
      a6 = adds[3];
 
    int a1o, a2o, a3o, a4o;
    find_linear_pairs(a1o,a2o,a3o,a4o,a1,a2,a3,a4,wtm,ic1);

//    ntor[0] = ntor[1] = 0;
//    nangles[0] = nangles[1] = 3;

    int nf = 0;
    int nang = 0;
   //add top/side: a1 shifted
    if (active1[a1])
    {
      int a2c = a2;
      if (a1o==a2)
        a2c = a3;

      ntor[nf+0] = ntor[nf+1] = 0;
      nangles[nf+0] = nangles[nf+1] = 4;

      angles[nang+0][0] = a5;
      angles[nang+0][1] = wtm;
      angles[nang+0][2] = a6;
      angles[nang+1][0] = a5;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a1o;
      angles[nang+2][0] = a1;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a6;
      angles[nang+3][0] = a1;
      angles[nang+3][1] = wtm;
      angles[nang+3][2] = a1o;

      anglev[nang+0] = 90.;
      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;
      anglev[nang+3] = 90.;

      angles[nang+4][0] = a5;
      angles[nang+4][1] = wtm;
      angles[nang+4][2] = a6;
      angles[nang+5][0] = a6;
      angles[nang+5][1] = wtm;
      angles[nang+5][2] = a1o;
      angles[nang+6][0] = a1;
      angles[nang+6][1] = wtm;
      angles[nang+6][2] = a5;
      angles[nang+7][0] = a1;
      angles[nang+7][1] = wtm;
      angles[nang+7][2] = a1o;

      anglev[nang+4] = 90.;
      anglev[nang+5] = 180.;
      anglev[nang+6] = 180.;
      anglev[nang+7] = 90.;

      nfound += 2;
      napp[nf+0] = napp[nf+1] = 2;
      get_vec_45p(a1,wtm,a2c,ic1.coords,&aprv[(1+2*nf)*3]);
      get_vec_45m(a1,wtm,a2c,ic1.coords,&aprv[(2+2*nf)*3]);
      for (int j=0;j<6;j++)
        aprv[(3+2*nf)*3+j] = aprv[(1+2*nf)*3+j];

      nf += 2;
      nang += 8;
    } //approach via a1
   //a2 shifted
    if (active1[a2])
    {
      int a2c = a3;
      if (a2o==a3)
        a2c = a1;

      ntor[nf+0] = ntor[nf+1] = 0;
      nangles[nf+0] = nangles[nf+1] = 4;

      angles[nang+0][0] = a5;
      angles[nang+0][1] = wtm;
      angles[nang+0][2] = a6;
      angles[nang+1][0] = a5;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a2o;
      angles[nang+2][0] = a2;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a6;
      angles[nang+3][0] = a2;
      angles[nang+3][1] = wtm;
      angles[nang+3][2] = a2o;

      anglev[nang+0] = 90.;
      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;
      anglev[nang+3] = 90.;

      angles[nang+4][0] = a5;
      angles[nang+4][1] = wtm;
      angles[nang+4][2] = a6;
      angles[nang+5][0] = a6;
      angles[nang+5][1] = wtm;
      angles[nang+5][2] = a2o;
      angles[nang+6][0] = a2;
      angles[nang+6][1] = wtm;
      angles[nang+6][2] = a5;
      angles[nang+7][0] = a2;
      angles[nang+7][1] = wtm;
      angles[nang+7][2] = a2o;

      anglev[nang+4] = 90.;
      anglev[nang+5] = 180.;
      anglev[nang+6] = 180.;
      anglev[nang+7] = 90.;

      nfound += 2;
      napp[nf+0] = napp[nf+1] = 2;
      get_vec_45p(a2,wtm,a2c,ic1.coords,&aprv[(1+2*nf)*3]);
      get_vec_45m(a2,wtm,a2c,ic1.coords,&aprv[(2+2*nf)*3]);
      for (int j=0;j<6;j++)
        aprv[(3+2*nf)*3+j] = aprv[(1+2*nf)*3+j];

      nf += 2;
      nang += 8;
    } //approach via a2
    if (active1[a3])
    {
      int a2c = a2;
      if (a3o==a2)
        a2c = a1;

      ntor[nf+0] = ntor[nf+1] = 0;
      nangles[nf+0] = nangles[nf+1] = 4;

      angles[nang+0][0] = a5;
      angles[nang+0][1] = wtm;
      angles[nang+0][2] = a6;
      angles[nang+1][0] = a5;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a3o;
      angles[nang+2][0] = a3;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a6;
      angles[nang+3][0] = a3;
      angles[nang+3][1] = wtm;
      angles[nang+3][2] = a3o;

      anglev[nang+0] = 90.;
      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;
      anglev[nang+3] = 90.;

      angles[nang+4][0] = a5;
      angles[nang+4][1] = wtm;
      angles[nang+4][2] = a6;
      angles[nang+5][0] = a6;
      angles[nang+5][1] = wtm;
      angles[nang+5][2] = a3o;
      angles[nang+6][0] = a3;
      angles[nang+6][1] = wtm;
      angles[nang+6][2] = a5;
      angles[nang+7][0] = a3;
      angles[nang+7][1] = wtm;
      angles[nang+7][2] = a3o;

      anglev[nang+4] = 90.;
      anglev[nang+5] = 180.;
      anglev[nang+6] = 180.;
      anglev[nang+7] = 90.;

      nfound += 2;
      napp[nf+0] = napp[nf+1] = 2;
      get_vec_45p(a3,wtm,a2c,ic1.coords,&aprv[(1+2*nf)*3]);
      get_vec_45m(a3,wtm,a2c,ic1.coords,&aprv[(2+2*nf)*3]);
      for (int j=0;j<6;j++)
        aprv[(3+2*nf)*3+j] = aprv[(1+2*nf)*3+j];

      nf += 2;
      nang += 8;
    } //approach via a3
    if (active1[a4])
    {
      int a2c = a2;
      if (a4o==a2)
        a2c = a3;

      ntor[nf+0] = ntor[nf+1] = 0;
      nangles[nf+0] = nangles[nf+1] = 4;

      angles[nang+0][0] = a5;
      angles[nang+0][1] = wtm;
      angles[nang+0][2] = a6;
      angles[nang+1][0] = a5;
      angles[nang+1][1] = wtm;
      angles[nang+1][2] = a4o;
      angles[nang+2][0] = a4;
      angles[nang+2][1] = wtm;
      angles[nang+2][2] = a6;
      angles[nang+3][0] = a4;
      angles[nang+3][1] = wtm;
      angles[nang+3][2] = a4o;

      anglev[nang+0] = 90.;
      anglev[nang+1] = 180.;
      anglev[nang+2] = 180.;
      anglev[nang+3] = 90.;

      angles[nang+4][0] = a5;
      angles[nang+4][1] = wtm;
      angles[nang+4][2] = a6;
      angles[nang+5][0] = a6;
      angles[nang+5][1] = wtm;
      angles[nang+5][2] = a4o;
      angles[nang+6][0] = a4;
      angles[nang+6][1] = wtm;
      angles[nang+6][2] = a5;
      angles[nang+7][0] = a4;
      angles[nang+7][1] = wtm;
      angles[nang+7][2] = a4o;

      anglev[nang+4] = 90.;
      anglev[nang+5] = 180.;
      anglev[nang+6] = 180.;
      anglev[nang+7] = 90.;

      nfound += 2;
      napp[nf+0] = napp[nf+1] = 2;
      get_vec_45p(a4,wtm,a2c,ic1.coords,&aprv[(1+2*nf)*3]);
      get_vec_45m(a4,wtm,a2c,ic1.coords,&aprv[(2+2*nf)*3]);
      for (int j=0;j<6;j++)
        aprv[(3+2*nf)*3+j] = aprv[(1+2*nf)*3+j];

      nf += 2;
      nang += 8;
    } //approach via a4

  } //nadd==2, nbrks==0
  else if (nadd==2 && nbrks==1)
  {
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1])
      a3 = btm[2][2];
    if (btm[3][2]==brks[0] || btm[3][2]==brks[1])
      a4 = btm[3][2];
    int a5 = adds[0];
    if (adds[0]==wtm)
      a5 = adds[1];
    int a6 = adds[2];
    if (adds[2]==wtm)
      a6 = adds[3];

    int a1l = brks[0]; //leaving atom
    if (brks[0]==wtm)
      a1l = brks[1];

    int a1o, a2o, a3o, a4o;
    find_linear_pairs(a1o,a2o,a3o,a4o,a1,a2,a3,a4,wtm,ic1);

    int a1lo = a1o; //opposite of leaving atom
    if (a1l==a2)
      a1lo = a2o;
    else if (a1l==a3)
      a1lo = a3o;
    else if (a1l==a4)
      a1lo = a4o;

    int a2c = a1; //perpendicular atom
    if (a1l!=a2 && a1lo!=a2)
      a2c = a2;
    if (a1l!=a3 && a1lo!=a3)
      a2c = a3;
    if (a1l!=a4 && a1lo!=a4)
      a2c = a4;

    ntor[0] = ntor[1] = 0;
    nangles[0] = nangles[1] = 2;

   //add one to top, one swap
    angles[0][0] = a5;
    angles[0][1] = wtm;
    angles[0][2] = a6;
    angles[1][0] = a5;
    angles[1][1] = wtm;
    angles[1][2] = a1lo;

    anglev[0] = 90.;
    anglev[1] = 180.;

    angles[2][0] = a5;
    angles[2][1] = wtm;
    angles[2][2] = a6;
    angles[3][0] = a6;
    angles[3][1] = wtm;
    angles[3][2] = a1lo;

    anglev[2] = 90.;
    anglev[3] = 180.;

    nfound = 2;
    napp[0] = napp[1] = 2;
    get_vec_45p(a1l,wtm,a2c,ic1.coords,&aprv[3]);
    get_vec_45m(a1l,wtm,a2c,ic1.coords,&aprv[6]);
    for (int j=0;j<6;j++)
      aprv[9+j] = aprv[3+j];
  } //nadd==2, nbrks==1
  else if (nadd==1 && nbrks==2)
  {
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1] || btm[0][2]==brks[2] || btm[0][2]==brks[3])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1] || btm[1][2]==brks[2] || btm[1][2]==brks[3])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1] || btm[2][2]==brks[2] || btm[2][2]==brks[3])
      a3 = btm[2][2];
    if (btm[3][2]==brks[0] || btm[3][2]==brks[1] || btm[3][2]==brks[2] || btm[3][2]==brks[3])
      a4 = btm[3][2];
    int a5 = adds[0];
    if (adds[0]==wtm)
      a5 = adds[1];

    int a1l = a3; //leaving group
    int a2l = a4; //leaving group
    int a1c = a1; //remaining group
    int a2c = a2; //remaining group 2
    if (brks[0]==a1 || brks[1]==a1 || brks[2]==a1 || brks[3]==a1)
    {
      a1l = a1;
      if (brks[0]==a2 || brks[1]==a2 || brks[2]==a2 || brks[3]==a2)
      {
        a2l = a2;
        a1c = a3;
        a2c = a4;
      }
      else
      {
        a1c = a2;
        if (brks[0]==a3 || brks[1]==a3 || brks[2]==a3 || brks[3]==a3)
        {
          a2l = a3;
          a2c = a4;
        }
        else
        {
          a2l = a4;
          a2c = a3;
        }
      }
    }
    else if (brks[0]==a2 || brks[1]==a2 || brks[2]==a2 || brks[3]==a2)
    {
      a1l = a2;
      a1c = a1;
      if (brks[0]==a3 || brks[1]==a3 || brks[2]==a3 || brks[3]==a3)
      {
        a2l = a3;
        a2c = a4;
      }
      else
      {
        a2l = a4;
        a2c = a3;
      }
    }
    else
    {
      a1c = a1;
      a2c = a2;
    }
    printf("   unbroken: %i %i \n",a1c+1,a2c+1);
    //printf("   a1-a5: %i %i %i %i %i \n",a1+1,a2+1,a3+1,a4+1,a5+1);
    //printf("   a1l: %i a2l: %i \n",a1l+1,a2l+1);

    ntor[0] = ntor[1] = 0;
    nangles[0] = nangles[1] = 1;

   //line up opposite
    angles[0][0] = a5;
    angles[0][1] = wtm;
    angles[0][2] = a1c;

    anglev[0] = 180.;

    angles[1][0] = a5;
    angles[1][1] = wtm;
    angles[1][2] = a2c;

    anglev[1] = 180.;

    nfound = 2;
    napp[0] = napp[1] = 2;
    get_vec_45p(a1l,wtm,a2c,ic1.coords,&aprv[3]);
    get_vec_45m(a1l,wtm,a2c,ic1.coords,&aprv[6]);
    get_vec_45p(a2l,wtm,a1c,ic1.coords,&aprv[9]);
    get_vec_45m(a2l,wtm,a1c,ic1.coords,&aprv[12]);
  }
  else
    nfound = 1;


  return nfound;
}


int ZStruct::run_3c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv)
{
  int nfound = 0;

  printf("  in run_3c_1 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  if (eta2_double_break(wtm,stm,3,nbrks,brks,ic1))
    nbrks--;

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;

  if (nadd==1 && nbrks==0) 
  {
    int a1 = btm[0][1];
    int a2 = btm[1][1];
    int a3 = btm[2][1];
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];

    double ang1 = ic1.angle_val(a1,wtm,a2);
    double ang2 = ic1.angle_val(a1,wtm,a3);
    double ang3 = ic1.angle_val(a2,wtm,a3);
    int a1c = a1;
    int a2c = a2;
    int a3c = a3;
    if (ang1>ang2 && ang1>ang3)
    {
      a1c = a3; a2c = a1; a3c = a2;
    }
    else if (ang2>ang1 && ang2>ang3)
    {
      a1c = a2; a2c = a1; a3c = a3;
    }
    else if (ang3>ang1 && ang3>ang2)
    {
      a1c = a1; a2c = a2; a3c = a3;
    }

    ntor[0] = 0;
    nangles[0] = 1;
    nangles[1] = nangles[2] = 2;

   //straight approach
    angles[0][0] = a1c;
    angles[0][1] = wtm;
    angles[0][2] = a4;

    anglev[0] = 180.;

   //bend in approach
    angles[1][0] = a2c;
    angles[1][1] = wtm;
    angles[1][2] = a4;
    angles[2][0] = a1c;
    angles[2][1] = wtm;
    angles[2][2] = a3c;

    anglev[1] = 180.;
    anglev[2] = 180.;

    angles[3][0] = a3c;
    angles[3][1] = wtm;
    angles[3][2] = a4;
    angles[4][0] = a1c;
    angles[4][1] = wtm;
    angles[4][2] = a2c;

    anglev[3] = 180.;
    anglev[4] = 180.;

    nfound = 3;

    napp[0] = napp[1] = napp[2] = 1;
   //approach opposite of a1c
    get_vec(a1c,wtm,ic1.coords,&aprv[3]);
   //Y involving a1c/aXc/wtm  X=2,3
    get_Y_vec(a1c,a3c,wtm,ic1.coords,&aprv[6]);
    get_Y_vec(a1c,a2c,wtm,ic1.coords,&aprv[9]);
    for (int j=6;j<11;j++) aprv[j] *= -1.;
  }
  else if (nadd==1 && nbrks==1) 
  {
    int a1 = btm[0][1];
    int a2 = btm[1][1];
    int a3 = btm[2][1];
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1])
      a3 = btm[2][2];
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];

    double ang1 = ic1.angle_val(a1,wtm,a2);
    double ang2 = ic1.angle_val(a1,wtm,a3);
    double ang3 = ic1.angle_val(a2,wtm,a3);
    int a1c = a1;
    int a2c = a2;
    int a3c = a3;
    if (ang1>ang2 && ang1>ang3)
    {
      a1c = a3; a2c = a1; a3c = a2;
    }
    else if (ang2>ang1 && ang2>ang3)
    {
      a1c = a2; a2c = a1; a3c = a3;
    }
    else if (ang3>ang1 && ang3>ang2)
    {
      a1c = a1; a2c = a2; a3c = a3;
    }
    int a2cc = a2c;
    int a3l = a3c;
    if (brks[0]==a2c || brks[1]==a2c)
    {
      a2cc = a3c;
      a3l = a2c;
    }

    ntor[0] = 0;
    nangles[0] = 1;
    nangles[1] = nangles[2] = 1;

   //straight approach
    angles[0][0] = a1c;
    angles[0][1] = wtm;
    angles[0][2] = a4;

    anglev[0] = 180.;

   //from above breaking group
    angles[1][0] = a2cc;
    angles[1][1] = wtm;
    angles[1][2] = a4;

    anglev[1] = 180.;

    nfound = 2;

    napp[0] = 1;
    napp[1] = 2;
    get_vec(a1c,wtm,ic1.coords,&aprv[3]);
    if (a1c==brks[0] || a1c==brks[1])
    {
      get_vec_45p(a3l,wtm,a2cc,ic1.coords,&aprv[6]);
      get_vec_45m(a3l,wtm,a2cc,ic1.coords,&aprv[9]);
    }
    else
    {
      get_vec_45p(a3l,wtm,a1c,ic1.coords,&aprv[6]);
      get_vec_45m(a3l,wtm,a1c,ic1.coords,&aprv[9]);
    }
  }
  else if (nadd==2 && nbrks==0) 
  {
    int a1 = btm[0][1];
    int a2 = btm[1][1];
    int a3 = btm[2][1];
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];
    int a5 = adds[2];
    if (adds[2]==wtm)
      a5 = adds[3];

    double ang1 = ic1.angle_val(a1,wtm,a2);
    double ang2 = ic1.angle_val(a1,wtm,a3);
    double ang3 = ic1.angle_val(a2,wtm,a3);
    int a1c = a1;
    int a2c = a2;
    int a3c = a3;
    if (ang1>ang2 && ang1>ang3)
    {
      a1c = a3; a2c = a1; a3c = a2;
    }
    else if (ang2>ang1 && ang2>ang3)
    {
      a1c = a2; a2c = a1; a3c = a3;
    }
    else if (ang3>ang1 && ang3>ang2)
    {
      a1c = a1; a2c = a2; a3c = a3;
    }

    ntor[0] = ntor[1] = ntor[2] = 0;
    nangles[0] = nangles[1] = nangles[2] = 3;

   //top and bottom add
    angles[0][0] = a4;
    angles[0][1] = wtm;
    angles[0][2] = a5;
    angles[1][0] = a1c;
    angles[1][1] = wtm;
    angles[1][2] = a4;
    angles[2][0] = a1c;
    angles[2][1] = wtm;
    angles[2][2] = a5;

    anglev[0] = 180.;
    anglev[1] = 90.;
    anglev[2] = 90.;

   //top and side add
    angles[3][0] = a4;
    angles[3][1] = wtm;
    angles[3][2] = a5;
    angles[4][0] = a1c;
    angles[4][1] = wtm;
    angles[4][2] = a4;
    angles[5][0] = a1c;
    angles[5][1] = wtm;
    angles[5][2] = a5;

    anglev[3] = 90.;
    anglev[4] = 180.;
    anglev[5] = 90.;

    angles[6][0] = a4;
    angles[6][1] = wtm;
    angles[6][2] = a5;
    angles[7][0] = a1c;
    angles[7][1] = wtm;
    angles[7][2] = a4;
    angles[8][0] = a1c;
    angles[8][1] = wtm;
    angles[8][2] = a5;

    anglev[6] = 90.;
    anglev[7] = 90.;
    anglev[8] = 180.;

    nfound = 3;

    napp[0] = napp[1] = napp[2] = 1;
    get_vec(a1c,wtm,ic1.coords,&aprv[3]);
    aprv[6] = aprv[3]; aprv[7] = aprv[4]; aprv[8] = aprv[5];
    aprv[9] = aprv[3]; aprv[10] = aprv[4]; aprv[11] = aprv[5];
  }
  else if (nadd==2 && nbrks==1) 
  {
    int a1 = btm[0][1];
    int a2 = btm[1][1];
    int a3 = btm[2][1];
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1])
      a3 = btm[2][2];
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];
    int a5 = adds[2];
    if (adds[2]==wtm)
      a5 = adds[3];

    double ang1 = ic1.angle_val(a1,wtm,a2);
    double ang2 = ic1.angle_val(a1,wtm,a3);
    double ang3 = ic1.angle_val(a2,wtm,a3);
    int a1c = a1;
    int a2c = a2;
    int a3c = a3;
    if (ang1>ang2 && ang1>ang3)
    {
      a1c = a3; a2c = a1; a3c = a2;
    }
    else if (ang2>ang1 && ang2>ang3)
    {
      a1c = a2; a2c = a1; a3c = a3;
    }
    else if (ang3>ang1 && ang3>ang2)
    {
      a1c = a1; a2c = a2; a3c = a3;
    }
    int a2cc = a2c;
    if (brks[0]==a2c || brks[1]==a2c)
      a2cc = a3c;

    ntor[0] = 0;
    nangles[0] = nangles[1] = 2;

   //straight approach
    angles[0][0] = a1c;
    angles[0][1] = wtm;
    angles[0][2] = a4;
    angles[1][0] = a2cc;
    angles[1][1] = wtm;
    angles[1][2] = a5;

    anglev[0] = 180.;
    anglev[1] = 180.;

    angles[2][0] = a1c;
    angles[2][1] = wtm;
    angles[2][2] = a5;
    angles[3][0] = a2cc;
    angles[3][1] = wtm;
    angles[3][2] = a4;

    anglev[2] = 180.;
    anglev[3] = 180.;

    nfound = 2;

    napp[0] = napp[1] = 1;
   //Y involving a1c/a2cc/wtm
    get_Y_vec(a1c,a2cc,wtm,ic1.coords,&aprv[3]);
    aprv[6] = aprv[3]; aprv[7] = aprv[4]; aprv[8] = aprv[5];
#if 0
    if (a1c==brks[0] || a1c==brks[1]) //handle second approach when a1c leaves
    {
      get_Y_vec(a2c,a3c,wtm,ic1.coords,&aprv[9]);
      aprv[12] = aprv[9]; aprv[13] = aprv[10]; aprv[14] = aprv[11];
      napp[0] = napp[1] = 2;
    }
#endif
  }
  else if (nadd==1 && nbrks==2) 
  {
    int a1 = btm[0][1];
    int a2 = btm[1][1];
    int a3 = btm[2][1];
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1] || btm[0][2]==brks[2] || btm[0][2]==brks[3])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1] || btm[1][2]==brks[2] || btm[1][2]==brks[3])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1] || btm[2][2]==brks[2] || btm[2][2]==brks[3])
      a3 = btm[2][2];
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];

    double ang1 = ic1.angle_val(a1,wtm,a2);
    double ang2 = ic1.angle_val(a1,wtm,a3);
    double ang3 = ic1.angle_val(a2,wtm,a3);
    int a1c = a1;
    int a2c = a2;
    int a3c = a3;
    if (ang1>ang2 && ang1>ang3)
    {
      a1c = a3; a2c = a1; a3c = a2;
    }
    else if (ang2>ang1 && ang2>ang3)
    {
      a1c = a2; a2c = a1; a3c = a3;
    }
    else if (ang3>ang1 && ang3>ang2)
    {
      a1c = a1; a2c = a2; a3c = a3;
    }
    int a0cc = a1c; //remaining group
    if (brks[0]==a1c || brks[1]==a1c || brks[2]==a1c || brks[3]==a1c)
    {
      a0cc = a2c;
      if (brks[0]==a2c || brks[1]==a2c || brks[2]==a2c || brks[3]==a2c)
        a0cc = a3c;
    }

    ntor[0] = 0;
    nangles[0] = nangles[1] = 1;

   //straight approach
    angles[0][0] = a0cc;
    angles[0][1] = wtm;
    angles[0][2] = a4;

    anglev[0] = 180.;

   //right angle approach
    angles[1][0] = a0cc;
    angles[1][1] = wtm;
    angles[1][2] = a4;

    anglev[1] = 90.;

    nfound = 2;

    napp[0] = napp[1] = 1;
    get_vec(a0cc,wtm,ic1.coords,&aprv[3]);
    aprv[6] = aprv[3]; aprv[7] = aprv[4]; aprv[8] = aprv[5];
  }
  else
    nfound = 1;



  return nfound;
}

int ZStruct::eta2_double_break(int wtm, int stm, int ncoord, int nbrks, int* brks, ICoord ic1)
{
  int found = 0;
  int b1,b2;

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;

  //if two breaks on an eta2 ligand, count as one break
  if (nbrks==2)
  {
    b1 = brks[0];
    b2 = brks[2];
    if (brks[0]==wtm)
      b1 = brks[1];
    if (brks[2]==wtm)
      b1 = brks[3];

    for (int i=0;i<ncoord;i++)
    if ((btm[i][1]==b1 && btm[i][2]==b2)
      || (btm[i][1]==b2 && btm[i][2]==b1))
    { 
      found = 1;
      break;
    }
  }

  if (found)
    printf("   b1: %i b2: %i are eta2 double break \n",b1+1,b2+1);

  return found;
}

int ZStruct::run_3c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv)
{
  int nfound = 0;

  printf("  in run_3c_0 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

  int** btm = ic1.bondstm;
  if (stm==2) 
    btm = ic1.bondstm2;
  int a1 = btm[0][1];
  int a2 = btm[1][1];
  int a3 = btm[2][1];

  if (eta2_double_break(wtm,stm,3,nbrks,brks,ic1))
    nbrks--;

  if (nadd==1 && nbrks==0) 
  {
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];

    ntor = 0;
    nangles[0] = nangles[1] = nangles[2] = nangles[3] = 1;

    angles[0][0] = a1;
    angles[0][1] = wtm;
    angles[0][2] = a4;
    angles[1][0] = a2;
    angles[1][1] = wtm;
    angles[1][2] = a4;
    angles[2][0] = a3;
    angles[2][1] = wtm;
    angles[2][2] = a4;
    angles[3][0] = a1;
    angles[3][1] = wtm;
    angles[3][2] = a4;

    anglev[0] = 160.;
    anglev[1] = 160.;
    anglev[2] = 160.;
    anglev[3] = 90.;

    nfound = 4;
    napp[0] = napp[1] = napp[2] = napp[3] = 1;
    get_vec(wtm,a1,ic1.coords,&aprv[3]);
    get_vec(wtm,a2,ic1.coords,&aprv[6]);
    get_vec(wtm,a3,ic1.coords,&aprv[9]);
   //note: only did 1 oop vec, could do 2
    get_oop_vec(a1,wtm,a2,ic1.coords,&aprv[12]);
  }
  else if (nadd==1 && nbrks==1) 
  {
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1])
      a3 = btm[2][2];
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];

    ntor = 0;
    nangles[0] = 2;

    int a1c = a1;
    int a2c = a2;
    int a3l = a3;
    if (brks[0]==a1 || brks[1]==a1)
    {
      a1c = a2;
      a2c = a3;
      a3l = a1;
    }
    else if (brks[0]==a2 || brks[1]==a2)
    {
      a1c = a1;
      a2c = a3;
      a3l = a2;
    }

    angles[0][0] = a1c;
    angles[0][1] = wtm;
    angles[0][2] = a4;
    angles[1][0] = a2c;
    angles[1][1] = wtm;
    angles[1][2] = a4;

    anglev[0] = 120.;
    anglev[1] = 120.;

    nfound = 1;
    napp[0] = 2;
    get_vec_45p(a3l,wtm,a1c,ic1.coords,&aprv[3]);
    get_vec_45m(a3l,wtm,a1c,ic1.coords,&aprv[6]);
  }
  else if (nadd==2 && nbrks==0) 
  {
    int a4 = adds[0];
    int a5 = adds[2];
    if (adds[0]==wtm)
      a4 = adds[1];
    if (adds[2]==wtm)
      a5 = adds[3];

    ntor[0] = 0;
    nangles[0] = 1;
 
   //ox add to top/bottom
    angles[0][0] = a4;
    angles[0][1] = wtm;
    angles[0][2] = a5;

    anglev[0] = 180.;

    ntor[1] = ntor[2] = ntor[3] = 0;
    nangles[1] = nangles[2] = nangles[3] = 2;

   //ox add to top/equitorial
    angles[1][0] = a4;
    angles[1][1] = wtm;
    angles[1][2] = a5;
    angles[2][0] = a1;
    angles[2][1] = wtm;
    angles[2][2] = a4;

    angles[3][0] = a4;
    angles[3][1] = wtm;
    angles[3][2] = a5;
    angles[4][0] = a2;
    angles[4][1] = wtm;
    angles[4][2] = a4;

    angles[5][0] = a4;
    angles[5][1] = wtm;
    angles[5][2] = a5;
    angles[6][0] = a3;
    angles[6][1] = wtm;
    angles[6][2] = a4;

    anglev[1] = 90.;
    anglev[2] = 180.;
    anglev[3] = 90.;
    anglev[4] = 180.;
    anglev[5] = 90.;
    anglev[6] = 180.;

   //note: doesn't distinguish "above" from "below"
    nfound = 4;
    napp[0] = napp[1] = napp[2] = napp[3] = 3;
   //3 approaches for type 1
    get_vec(wtm,a1,ic1.coords,&aprv[3]);
    get_vec(wtm,a2,ic1.coords,&aprv[6]);
    get_vec(wtm,a3,ic1.coords,&aprv[9]);
   //3 approaches for type 2 (same as above)
    for (int j=0;j<27;j++)
      aprv[12+j] = aprv[3+j];
  }
  else if (nadd==2 && nbrks==1) 
  {
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1])
      a3 = btm[2][2];
    int a4 = adds[0];
    int a5 = adds[2];
    if (adds[0]==wtm)
      a4 = adds[1];
    if (adds[2]==wtm)
      a5 = adds[3];

    int a1c = a1;
    int a2c = a2;
    if (brks[0]==a1 || brks[1]==a1)
    {
      a1c = a2;
      a2c = a3;
    }
    else if (brks[0]==a2 || brks[1]==a2)
    {
      a1c = a1;
      a2c = a3;
    }

    ntor[0] = ntor[1] = 0;
    nangles[0] = nangles[1] = 2;

    angles[0][0] = a1c;
    angles[0][1] = wtm;
    angles[0][2] = a4;
    angles[1][0] = a2c;
    angles[1][1] = wtm;
    angles[1][2] = a5;

    angles[2][0] = a1c;
    angles[2][1] = wtm;
    angles[2][2] = a5;
    angles[3][0] = a2c;
    angles[3][1] = wtm;
    angles[3][2] = a4;

    anglev[0] = 100.;
    anglev[1] = 100.;
    anglev[2] = 100.;
    anglev[3] = 100.;

    nfound = 2;
    napp[0] = napp[1] = 4;
    get_vec(wtm,a1c,ic1.coords,&aprv[3]);
    get_vec(wtm,a2c,ic1.coords,&aprv[6]);
    get_Y_vec_45p(a1c,a2c,wtm,ic1.coords,&aprv[9]);
    get_Y_vec_45m(a1c,a2c,wtm,ic1.coords,&aprv[12]);
    for (int j=0;j<12;j++)
      aprv[15+j] = aprv[3+j];
  }
  else if (nadd==1 && nbrks==2) 
  {
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1] || btm[0][2]==brks[2] || btm[0][2]==brks[3])
      a1 = btm[0][2];
    if (btm[1][2]==brks[0] || btm[1][2]==brks[1] || btm[1][2]==brks[2] || btm[1][2]==brks[3])
      a2 = btm[1][2];
    if (btm[2][2]==brks[0] || btm[2][2]==brks[1] || btm[2][2]==brks[2] || btm[2][2]==brks[3])
      a3 = btm[2][2];
    int a4 = adds[0];
    if (adds[0]==wtm)
      a4 = adds[1];

    // note: doesn't handle 2 breaks on one eta2 (but should never get here)
    int a1c = a1; //remaining group
    int a2l = a2;
    int a3l = a3;
    if (brks[0]==a1 || brks[1]==a1 || brks[2]==a1 || brks[3]==a1)
    {
      a1c = a2;
      a2l = a1;
      a3l = a3;
      if (brks[0]==a2 || brks[1]==a2 || brks[2]==a2 || brks[3]==a2)
      {
        a1c = a3;
        a2l = a1;
        a3l = a2;
      }
    }

    ntor[0] = 0;
    nangles[0] = 1;

    angles[0][0] = a1c;
    angles[0][1] = wtm;
    angles[0][2] = a4;

    anglev[0] = 120.;

    nfound = 1;
    napp[0] = 2;
    get_vec(a2l,wtm,ic1.coords,&aprv[3]);
    get_vec(a3l,wtm,ic1.coords,&aprv[6]);
  }
  else
    nfound = 1;


  return nfound;
}

int ZStruct::run_2c_0(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv)
{
  int nfound = 0;

  printf("  in run_2c_0 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;
  if (nbrks==2 && nadd==1) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];
 
 //no nbrks==2 here

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;

  if (nadd>0) 
  {
    //find possible Y
    int a1 = btm[0][1];
    int a2 = btm[1][1];
    if (nbrks>0)
    {
      if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
        a1 = btm[0][2];
      if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
        a2 = btm[1][2];
    }
    int a3 = adds[0];
    if (adds[0]==wtm)
      a3 = adds[1];

    if (nadd==1 && nbrks==1)
    {
      int a1c = a1; //angle with staying group
      int a0l = a2;
      if (a1==brks[0] || a1==brks[1])
      {
        a1c = a2;
        a0l = a1;
      }

      ntor[0] = 0;
      nangles[0] = 1;
      angles[0][0] = a1c;
      angles[0][1] = wtm;
      angles[0][2] = a3;

      anglev[0] = 120.;

      nfound = 1;
      napp[0] = 3;
      get_Y_vec(a1,a2,wtm,ic1.coords,&aprv[3]);
      get_vec_45p(a0l,wtm,a1c,ic1.coords,&aprv[6]);
      get_vec_45m(a0l,wtm,a1c,ic1.coords,&aprv[9]);
    }
    else if (nadd==1 && nbrks==0)
    {
      ntor[0] = 0;
      nangles[0] = 2;
      angles[0][0] = a1;
      angles[0][1] = wtm;
      angles[0][2] = a3;
      angles[1][0] = a2;
      angles[1][1] = wtm;
      angles[1][2] = a3;

      anglev[0] = 120.; //create Y
      anglev[1] = 120.;

      nfound = 1;
      napp[0] = 1;
      get_Y_vec(a1,a2,wtm,ic1.coords,&aprv[3]);
    }
    else if (nadd==2 && nbrks==1)
    {
      int a1c = a1; //staying group
      if (a1==brks[0] || a1==brks[1])
        a1c = a2;
      int a4 = adds[2];
      if (adds[2]==wtm)
        a4 = adds[3];

      ntor[0] = 0;
      nangles[0] = 2;
      angles[0][0] = a1c;
      angles[0][1] = wtm;
      angles[0][2] = a3;
      angles[1][0] = a1c;
      angles[1][1] = wtm;
      angles[1][2] = a4;

      anglev[0] = 110.;
      anglev[1] = 110.;

      nfound = 1;
      napp[0] = 1; 
      get_vec(a1c,wtm,ic1.coords,&aprv[3]);
    }
    else if (nadd==2 && nbrks==0)
    {
      ntor[0] = 0;
      nangles[0] = 1;
      angles[0][0] = a1;
      angles[0][1] = wtm;
      angles[0][2] = a2;

      anglev[0] = 110.;

      nfound = 1;
      napp[0] = 1;
      get_Y_vec(a1,a2,wtm,ic1.coords,&aprv[3]);
    }
  }
  else
    nfound = 1;


  return nfound;
}

int ZStruct::run_2c_1(int wtm, int stm, int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int* nangles, int** angles, double* anglev, int* ntor, int** tor, double* torv, int* napp, double* aprv)
{
  int nfound = 0;

  printf("  in run_2c_1 for naddtm: %i nbrktm: %i \n",nadd,nbrks);

  ntor[0] = 0;
  nangles[0] = 0;
  napp[0] = 1;

 //not handling these cases
  if (nbrks>2 || nadd>2) return 0;
  if (nbrks==2 && nadd==1) return 0;

 //center align vector at metal
  aprv[0] = ic1.coords[3*wtm+0];
  aprv[1] = ic1.coords[3*wtm+1];
  aprv[2] = ic1.coords[3*wtm+2];

 //no nbrks==2 here

  int** btm = ic1.bondstm;
  if (stm==2)
    btm = ic1.bondstm2;

  if (nadd==1) //0 or 1 break
  {
    //find possible "T"
    int a1 = btm[0][1];
    int a2 = btm[1][1];
    if (nbrks>0)
    {
      if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
        a1 = btm[0][2];
      if (btm[1][2]==brks[0] || btm[1][2]==brks[1])
        a2 = btm[1][2];
    }
    if (nbrks==1)
    if (brks[0]==a1 || brks[1]==a1) //angle must be with remaining group
    {
      int tmp = a2;
      a2 = a1;
      a1 = tmp;
    }
    int a3 = adds[0];
    if (adds[0]==wtm)
      a3 = adds[1];

    ntor[0] = 0;
    nangles[0] = 1;
    angles[0][0] = a1;
    angles[0][1] = wtm;
    angles[0][2] = a3;

    printf("   nadd: 1. angles: %i %i %i \n",a1+1,wtm+1,a3+1);

    anglev[0] = 90.; //add to create T
    if (nbrks==1) anglev[0] = 180.; //swap move

    nfound = 1;
    napp[0] = 2; 
    //use cross product of a1-->wtm...
    get_perp_vec(a1,wtm,ic1.coords,&aprv[3]);
    aprv[6] = -aprv[3]; aprv[7] = -aprv[4]; aprv[8] = -aprv[5];
//    printf("   NOTE: need to do 3 rotations of this type \n");
#if 0
    printf("   finding perpendicular v to %i %i \n",a1,wtm);
    print_xyz_gen(ic1.natoms,ic1.anames,ic1.coords);
    printf("    X %4.3f %4.3f %4.3f \n",aprv[0],aprv[1],aprv[2]);
    printf("    X %4.3f %4.3f %4.3f \n",aprv[0]+aprv[3],aprv[1]+aprv[4],aprv[2]+aprv[5]);
#endif
#if 0
    if (nbrks==1) //Not using
    {
      napp[0] += 2;
      get_vec_45p(a1,wtm,a2,ic1.coords,&aprv[6]);
      get_vec_45m(a1,wtm,a2,ic1.coords,&aprv[9]);
    }
#endif
  }
  else if (nadd==2 && nbrks==0)
  {
    //find attached
    int a1 = btm[0][1];
    int a2 = btm[1][1];

    ntor[0] = 0;
    nangles[0] = 1;
    angles[0][0] = angles[2][0] = a1;
    angles[0][1] = angles[2][1] = wtm;
    angles[0][2] = angles[2][2] = a2;

    anglev[0] = 110.; //bend back to make room

    nfound = 1;
    napp[0] = 1;
    //use cross product of a1-->wtm...
    get_perp_vec(a1,wtm,ic1.coords,&aprv[3]);
    aprv[6] = -aprv[3]; aprv[7] = -aprv[4]; aprv[8] = -aprv[5];
//    printf("   NOTE: need to do 3 rotations of this type \n");
  }
  else if (nadd==2 && nbrks==1)
  {
    //find possible "T"
    int a1 = btm[0][1];
    if (btm[0][2]==brks[0] || btm[0][2]==brks[1])
      a1 = btm[0][2];
    if (brks[0]==a1 || brks[1]==a1) //angle must be with remaining group
      a1 = btm[1][1];
    int a2a = adds[0];
    if (adds[0]==wtm)
      a2a = adds[1];
    int a2b = adds[2];
    if (adds[2]==wtm)
      a2b = adds[3];

    ntor[0] = ntor[1] = 0;
    nangles[0] = nangles[1] = 2;
    angles[0][0] = angles[2][0] = a1;
    angles[0][1] = angles[2][1] = wtm;
    angles[0][2] = angles[2][2] = a2a;
    angles[1][0] = angles[3][0] = a1;
    angles[1][1] = angles[3][1] = wtm;
    angles[1][2] = angles[3][2] = a2b;

    anglev[0] = 90.; //creates T
    anglev[1] = 180.; 
    anglev[2] = 180.; //creates T, attaching swapped
    anglev[3] = 90.; 

    nfound = 2;
    napp[0] = napp[1] = 2;
    //use cross product of a1-->wtm...
    get_perp_vec(a1,wtm,ic1.coords,&aprv[3]);
    aprv[6] = -aprv[3]; aprv[7] = -aprv[4]; aprv[8] = -aprv[5];
    for (int j=0;j<6;j++)
      aprv[9+j] = aprv[3+j];
//    printf("   NOTE: need to do 3 rotations of this type \n");
  }
  else
    nfound = 1;


  return nfound;
}


#if 0
//doesn't work
void ZStruct::align_vec(int a1, int a2, double* coords, double* w)
{
  double x1 = coords[3*a1+0] - coords[3*a2+0];
  double x2 = coords[3*a1+1] - coords[3*a2+1];
  double x3 = coords[3*a1+2] - coords[3*a2+2];
  double n1 = sqrt(x1*x1+x2*x2+x3*x3);
  x1 = x1/n1; x2=x2/n1; x3=x3/n1;

  double dot = w[0]*x1 + w[1]*x2 + w[2]*x3;
  if (dot<0.)
  {
    w[0] *= -1;
    w[1] *= -1;
    w[2] *= -1;
  }

  return;
}
#endif

//a2 points to a1, 45 degrees above a1-a2-a3 plane (normal vector points to a4 from a2)
void ZStruct::get_vec_45p_align(int a1, int a2, int a3, int a4, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_vec_45p_align requires three unique atoms \n");
    exit(1);
  }

  //NOTE: could generalize for any angle

  w[0] = coords[3*a1+0] - coords[3*a2+0];
  w[1] = coords[3*a1+1] - coords[3*a2+1];
  w[2] = coords[3*a1+2] - coords[3*a2+2];

  double n1 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n1;
  w[1] = w[1]/n1;
  w[2] = w[2]/n1;

  double* u = new double[3];
  double* v = new double[3];
  double* z = new double[3];
  u[0] = u[1] = u[2] = 0.;
  v[0] = v[1] = v[2] = 0.;
  z[0] = z[1] = z[2] = 0.;

  u[0] = coords[3*a3+0] - coords[3*a2+0];
  u[1] = coords[3*a3+1] - coords[3*a2+1];
  u[2] = coords[3*a3+2] - coords[3*a2+2];

  double n2 = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  u[0] = u[0]/n2;
  u[1] = u[1]/n2;
  u[2] = u[2]/n2;

  cross (v,u,w); //v is resultant vector

  double n3 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n3;
  v[1] = v[1]/n3;
  v[2] = v[2]/n3;

  z[0] = coords[3*a4+0] - coords[3*a2+0];
  z[1] = coords[3*a4+1] - coords[3*a2+1];
  z[2] = coords[3*a4+2] - coords[3*a2+2];

  double dot = v[0]*z[0] + v[1]*z[1] + v[2]*z[2];
  if (dot<0.)
  {
    v[0] *= -1.;
    v[1] *= -1.;
    v[2] *= -1.;
  }

  //point w 45 degrees above plane
  w[0] += v[0];
  w[1] += v[1];
  w[2] += v[2];

  double n4 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n4;
  w[1] = w[1]/n4;
  w[2] = w[2]/n4;

  delete [] u;
  delete [] v;
  delete [] z;

  return;
}

//a2 points to a1, 45 degrees above a1-a2-a3 plane (normal vector points to a4 from a2)
void ZStruct::get_vec_45m_align(int a1, int a2, int a3, int a4, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_vec_45m requires three unique atoms \n");
    exit(1);
  }

  //NOTE: could generalize for any angle

  w[0] = coords[3*a1+0] - coords[3*a2+0];
  w[1] = coords[3*a1+1] - coords[3*a2+1];
  w[2] = coords[3*a1+2] - coords[3*a2+2];

  double n1 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n1;
  w[1] = w[1]/n1;
  w[2] = w[2]/n1;

  double* u = new double[3];
  double* v = new double[3];
  double* z = new double[3];
  u[0] = u[1] = u[2] = 0.;
  v[0] = v[1] = v[2] = 0.;
  z[0] = z[1] = z[2] = 0.;

  u[0] = coords[3*a3+0] - coords[3*a2+0];
  u[1] = coords[3*a3+1] - coords[3*a2+1];
  u[2] = coords[3*a3+2] - coords[3*a2+2];

  double n2 = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  u[0] = u[0]/n2;
  u[1] = u[1]/n2;
  u[2] = u[2]/n2;

  cross (v,u,w); //v is resultant vector

  double n3 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n3;
  v[1] = v[1]/n3;
  v[2] = v[2]/n3;

  z[0] = coords[3*a4+0] - coords[3*a2+0];
  z[1] = coords[3*a4+1] - coords[3*a2+1];
  z[2] = coords[3*a4+2] - coords[3*a2+2];

  double dot = v[0]*z[0] + v[1]*z[1] + v[2]*z[2];
  if (dot>0.)
  {
    v[0] *= -1.;
    v[1] *= -1.;
    v[2] *= -1.;
  }

  //point w 45 degrees above plane
  w[0] += v[0];
  w[1] += v[1];
  w[2] += v[2];

  double n4 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n4;
  w[1] = w[1]/n4;
  w[2] = w[2]/n4;

  delete [] u;
  delete [] v;
  delete [] z;

  return;
}


//a2 points to a1, 45 degrees above a1-a2-a3 plane
void ZStruct::get_vec_45p(int a1, int a2, int a3, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_vec_45p requires three unique atoms \n");
    exit(1);
  }

  //NOTE: could generalize for any angle

  w[0] = coords[3*a1+0] - coords[3*a2+0];
  w[1] = coords[3*a1+1] - coords[3*a2+1];
  w[2] = coords[3*a1+2] - coords[3*a2+2];

  double n1 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n1;
  w[1] = w[1]/n1;
  w[2] = w[2]/n1;

  double* u = new double[3];
  double* v = new double[3];
  u[0] = u[1] = u[2] = 0.;
  v[0] = v[1] = v[2] = 0.;

  u[0] = coords[3*a3+0] - coords[3*a2+0];
  u[1] = coords[3*a3+1] - coords[3*a2+1];
  u[2] = coords[3*a3+2] - coords[3*a2+2];

  double n2 = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  u[0] = u[0]/n2;
  u[1] = u[1]/n2;
  u[2] = u[2]/n2;

  cross (v,u,w); //v is resultant vector

  double n3 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n3;
  v[1] = v[1]/n3;
  v[2] = v[2]/n3;

  //point w 45 degrees above plane
  w[0] += v[0];
  w[1] += v[1];
  w[2] += v[2];

  double n4 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n4;
  w[1] = w[1]/n4;
  w[2] = w[2]/n4;

  delete [] u;
  delete [] v;

  return;
}

void ZStruct::get_vec_45m(int a1, int a2, int a3, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_vec_45m requires three unique atoms \n");
    exit(1);
  }

  w[0] = coords[3*a1+0] - coords[3*a2+0];
  w[1] = coords[3*a1+1] - coords[3*a2+1];
  w[2] = coords[3*a1+2] - coords[3*a2+2];

  double n1 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n1;
  w[1] = w[1]/n1;
  w[2] = w[2]/n1;

  double* u = new double[3];
  double* v = new double[3];
  u[0] = u[1] = u[2] = 0.;
  v[0] = v[1] = v[2] = 0.;

  u[0] = coords[3*a3+0] - coords[3*a2+0];
  u[1] = coords[3*a3+1] - coords[3*a2+1];
  u[2] = coords[3*a3+2] - coords[3*a2+2];

  double n2 = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  u[0] = u[0]/n2;
  u[1] = u[1]/n2;
  u[2] = u[2]/n2;

  cross (v,u,w); //v is resultant vector

  double n3 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n3;
  v[1] = v[1]/n3;
  v[2] = v[2]/n3;

  //point w 45 degrees below plane
  w[0] -= v[0];
  w[1] -= v[1];
  w[2] -= v[2];

  double n4 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n4;
  w[1] = w[1]/n4;
  w[2] = w[2]/n4;

  delete [] u;
  delete [] v;

  return;
}

//a2 points to a1
void ZStruct::get_vec(int a1, int a2, double* coords, double* w)
{
  if (a1==a2)
  {
    printf("\n ERROR: get_vec requires two unique atoms \n");
    exit(1);
  }

  w[0] = coords[3*a1+0] - coords[3*a2+0];
  w[1] = coords[3*a1+1] - coords[3*a2+1];
  w[2] = coords[3*a1+2] - coords[3*a2+2];

  double n1 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n1;
  w[1] = w[1]/n1;
  w[2] = w[2]/n1;

  return;
}

void ZStruct::get_Y_vec(int a1, int a2, int a3, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_Y_vec requires three unique atoms \n");
    exit(1);
  }

  double x1 = coords[3*a1+0] - coords[3*a3+0];
  double x2 = coords[3*a1+1] - coords[3*a3+1];
  double x3 = coords[3*a1+2] - coords[3*a3+2];
  double n1 = sqrt(x1*x1+x2*x2+x3*x3);
  x1 = x1/n1; x2=x2/n1; x3=x3/n1;

  double y1 = coords[3*a2+0] - coords[3*a3+0];
  double y2 = coords[3*a2+1] - coords[3*a3+1];
  double y3 = coords[3*a2+2] - coords[3*a3+2];
  double n2 = sqrt(y1*y1+y2*y2+y3*y3);
  y1 = y1/n2; y2=y2/n2; y3=y3/n2;

  w[0] = - (x1+y1)/2.;
  w[1] = - (x2+y2)/2.;
  w[2] = - (x3+y3)/2.;

  double n3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n3;
  w[1] = w[1]/n3;
  w[2] = w[2]/n3;

  return;
}

void ZStruct::get_Y_vec_45p(int a1, int a2, int a3, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_perp_vec requires three unique atoms \n");
    exit(1);
  }

  double* u = new double[3];
  double* v = new double[3];

  u[0] = coords[3*a1+0] - coords[3*a3+0];
  u[1] = coords[3*a1+1] - coords[3*a3+1];
  u[2] = coords[3*a1+2] - coords[3*a3+2];
  double n1 = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  u[0] = u[0]/n1; u[1]=u[1]/n1; u[2]=u[2]/n1;

  v[0] = coords[3*a2+0] - coords[3*a3+0];
  v[1] = coords[3*a2+1] - coords[3*a3+1];
  v[2] = coords[3*a2+2] - coords[3*a3+2];
  double n2 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n2; v[1]=v[1]/n2; v[2]=v[2]/n2;

  w[0] = - (u[0]+v[0])/2.;
  w[1] = - (u[1]+v[1])/2.;
  w[2] = - (u[2]+v[2])/2.;

  double n3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n3;
  w[1] = w[1]/n3;
  w[2] = w[2]/n3;

  cross (v,u,w); //v is resultant vector

  n3 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n3;
  v[1] = v[1]/n3;
  v[2] = v[2]/n3;

  //point w 45 degrees above plane
  w[0] += v[0];
  w[1] += v[1];
  w[2] += v[2];

  double n4 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n4;
  w[1] = w[1]/n4;
  w[2] = w[2]/n4;

  delete [] u;
  delete [] v;

  return;
}

void ZStruct::get_Y_vec_45m(int a1, int a2, int a3, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_perp_vec requires three unique atoms \n");
    exit(1);
  }

  double* u = new double[3];
  double* v = new double[3];

  u[0] = coords[3*a1+0] - coords[3*a3+0];
  u[1] = coords[3*a1+1] - coords[3*a3+1];
  u[2] = coords[3*a1+2] - coords[3*a3+2];
  double n1 = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  u[0] = u[0]/n1; u[1]=u[1]/n1; u[2]=u[2]/n1;

  v[0] = coords[3*a2+0] - coords[3*a3+0];
  v[1] = coords[3*a2+1] - coords[3*a3+1];
  v[2] = coords[3*a2+2] - coords[3*a3+2];
  double n2 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n2; v[1]=v[1]/n2; v[2]=v[2]/n2;

  w[0] = - (u[0]+v[0])/2.;
  w[1] = - (u[1]+v[1])/2.;
  w[2] = - (u[2]+v[2])/2.;

  double n3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n3;
  w[1] = w[1]/n3;
  w[2] = w[2]/n3;

  cross (v,u,w); //v is resultant vector

  n3 = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = v[0]/n3;
  v[1] = v[1]/n3;
  v[2] = v[2]/n3;

  //point w 45 degrees above plane
  w[0] -= v[0];
  w[1] -= v[1];
  w[2] -= v[2];

  double n4 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0]/n4;
  w[1] = w[1]/n4;
  w[2] = w[2]/n4;

  delete [] u;
  delete [] v;

  return;
}

void ZStruct::get_oop_vec(int a1, int a2, int a3, double* coords, double* w)
{
  if (a1==a2 || a1==a3 || a2==a3)
  {
    printf("\n ERROR: get_perp_vec requires three unique atoms \n");
    exit(1);
  }

  double* x1 = new double[3];
  x1[0] = coords[3*a1+0] - coords[3*a2+0];
  x1[1] = coords[3*a1+1] - coords[3*a2+1];
  x1[2] = coords[3*a1+2] - coords[3*a2+2];

  double n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
  x1[0] = x1[0] / n1;
  x1[1] = x1[1] / n1;
  x1[2] = x1[2] / n1;

  double* x2 = new double[3];
  x2[0] = coords[3*a3+0] - coords[3*a2+0];
  x2[1] = coords[3*a3+1] - coords[3*a2+1];
  x2[2] = coords[3*a3+2] - coords[3*a2+2];

  double n2 = sqrt(x2[0]*x2[0]+x2[1]*x2[1]+x2[2]*x2[2]);
  x2[0] = x2[0] / n2;
  x2[1] = x2[1] / n2;
  x2[2] = x2[2] / n2;

  cross(w,x1,x2);
  double n3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  w[0] = w[0] / n3;
  w[1] = w[1] / n3;
  w[2] = w[2] / n3;

  delete [] x1;
  delete [] x2;

  return;
}

void ZStruct::get_perp_vec(int a1, int a2, double* coords, double* w)
{
  if (a1==a2)
  {
    printf("\n ERROR: get_perp_vec requires two unique atoms \n");
    exit(1);
  }

  double* x1 = new double[3];
  x1[0] = coords[3*a1+0] - coords[3*a2+0];
  x1[1] = coords[3*a1+1] - coords[3*a2+1];
  x1[2] = coords[3*a1+2] - coords[3*a2+2];
  double* x2 = new double[3];
  x2[0] = x2[1] = 0.;
  x2[2] = 1.;

  double n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
  x1[0] = x1[0] / n1;
  x1[1] = x1[1] / n1;
  x1[2] = x1[2] / n1;

  cross (w,x1,x2);
 
  double n2 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  if (n2<0.00000000001)
  {
    printf("\n ERROR: cannot make perpendicular vector \n");
    printf(" values: %8.6f %8.6f %8.6f \n",w[0],w[1],w[2]);
    exit(1);
  }

  w[0] = w[0]/n2;
  w[1] = w[1]/n2;
  w[2] = w[2]/n2;


  delete [] x1;
  delete [] x2;

  return;
}


int ZStruct::find_2c(int wtm, int stm, ICoord ic1)
{
  int type = 0;

  int* attached = new int[4];
  int nattached = find_attached(wtm,stm,ic1,attached);

  double angle = 0.;
  if (nattached==2)
    angle = ic1.angle_val(attached[0],wtm,attached[1]);
  else if (nattached==3)
  {
    printf("   need to handle eta2 \n");
    int a1 = attached[2];
    int a2a = attached[0];
    int a2b = attached[1];

    angle = ic1.angle_val_eta2(a1,wtm,a2a,a2b);
  }
  else if (nattached==4)
  {
    printf("   need to handle double eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];

    angle = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
  }
  else if (nattached>4)
  {
    printf("\n  ERROR: greater than eta2 not implemented \n");
    exit(1);
  }

  printf("   angle: %5.2f \n",angle);

  if (angle>ANGLE_LINEAR) //linear
    type = 1;
  else            //bent
    type = 0;

  delete [] attached;

  return type;
}

int ZStruct::find_3c(int wtm, int stm, ICoord ic1)
{
  int type = 0;

  int* attached = new int[3];
  int nattached = find_attached(wtm,stm,ic1,attached);

  double* angles = new double[6];
  if (nattached==3)
  {
    angles[0] = ic1.angle_val(attached[0],wtm,attached[1]);
    angles[1] = ic1.angle_val(attached[0],wtm,attached[2]);
    angles[2] = ic1.angle_val(attached[1],wtm,attached[2]);
  }
  else if (nattached==4)
  {
    printf("   need to handle single eta2 \n");
    int a1 = attached[2];
    int a2 = attached[3];
    int a3a = attached[0];
    int a3b = attached[1];

    angles[0] = ic1.angle_val(a1,wtm,a2);
    angles[1] = ic1.angle_val_eta2(a1,wtm,a3a,a3b);
    angles[2] = ic1.angle_val_eta2(a2,wtm,a3a,a3b);
  }
  else if (nattached==5)
  {
    printf("   need to handle double eta2 \n");
    int a1 = attached[4];
    int a2a = attached[0];
    int a2b = attached[1];
    int a3a = attached[2];
    int a3b = attached[3];

    angles[0] = ic1.angle_val_eta2(a1,wtm,a2a,a2b);
    angles[1] = ic1.angle_val_eta2(a1,wtm,a3a,a3b);
    angles[2] = ic1.angle_val_eta2_eta2(a2a,a2b,wtm,a3a,a3b);
  }
  else if (nattached==6)
  {
    printf("   need to handle triple eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];
    int a3a = attached[4];
    int a3b = attached[5];

    angles[0] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
    angles[1] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a3a,a3b);
    angles[2] = ic1.angle_val_eta2_eta2(a2a,a2b,wtm,a3a,a3b);
  }
  else if (nattached>6)
  {
    printf("\n  ERROR: greater than eta2 not implemented \n");
    exit(1);
  }

  sort_angles(3,angles);
  printf("   angle1: %5.2f angle2: %5.2f angle3: %5.2f \n",angles[0],angles[1],angles[2]);

  if (angles[0]>ANGLE_LINEAR) //T shaped
    type = 1;
  else            //trigonal
    type = 0;

  delete [] angles;
  delete [] attached;

  return type;
}


int ZStruct::find_4c(int wtm, int stm, ICoord ic1)
{
  int type = 0;

  int* attached = new int[8];
  int nattached = find_attached(wtm,stm,ic1,attached);

  double* angles = new double[6];
  if (nattached==4)
  {
    angles[0] = ic1.angle_val(attached[0],wtm,attached[1]);
    angles[1] = ic1.angle_val(attached[0],wtm,attached[2]);
    angles[2] = ic1.angle_val(attached[0],wtm,attached[3]);
    angles[3] = ic1.angle_val(attached[1],wtm,attached[2]);
    angles[4] = ic1.angle_val(attached[1],wtm,attached[3]);
    angles[5] = ic1.angle_val(attached[2],wtm,attached[3]);
  }
  else if (nattached==5)
  {
    printf("   need to handle single eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2 = attached[2];
    int a3 = attached[3];
    int a4 = attached[4];

    angles[0] = ic1.angle_val(a2,wtm,a3);
    angles[1] = ic1.angle_val(a2,wtm,a4);
    angles[2] = ic1.angle_val(a3,wtm,a4);
    angles[3] = ic1.angle_val_eta2(a2,wtm,a1a,a1b);
    angles[4] = ic1.angle_val_eta2(a3,wtm,a1a,a1b);
    angles[5] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
  }
  else if (nattached==6)
  {
    printf("   need to handle double eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];
    int a3 = attached[4];
    int a4 = attached[5];

    angles[0] = ic1.angle_val(a3,wtm,a4);
    angles[1] = ic1.angle_val_eta2(a3,wtm,a1a,a1b);
    angles[2] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[3] = ic1.angle_val_eta2(a3,wtm,a2a,a2b);
    angles[4] = ic1.angle_val_eta2(a4,wtm,a2a,a2b);
    angles[5] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
  }
  else if (nattached==7)
  {
    printf("   need to handle triple eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];
    int a3a = attached[4];
    int a3b = attached[5];
    int a4 = attached[6];

    angles[0] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[1] = ic1.angle_val_eta2(a4,wtm,a2a,a2b);
    angles[2] = ic1.angle_val_eta2(a4,wtm,a3a,a3b);
    angles[3] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
    angles[4] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a3a,a3b);
    angles[5] = ic1.angle_val_eta2_eta2(a2a,a2b,wtm,a3a,a3b);
  }
  else 
  {
    printf("\n ERROR: not yet handling quadruple eta2 \n");
    exit(1);
  }

  sort_angles(6,angles);
  printf("   angle1: %5.2f angle2: %5.2f angle3: %5.2f angle4: %5.2f angle5: %5.2f angle6: %5.2f \n",angles[0],angles[1],angles[2],angles[3],angles[4],angles[5]);

  if (angles[0]>ANGLE_LINEAR && angles[1]>ANGLE_LINEAR) //square planar
    type = 1;
  else            //tetrahedral
    type = 0;

  delete [] angles;
  delete [] attached;

  return type;
}


int ZStruct::find_5c(int wtm, int stm, ICoord ic1)
{
  int type = 0;

  int* attached = new int[10];
  int nattached = find_attached(wtm,stm,ic1,attached);

  double* angles = new double[10];
  find_5c_angles(nattached,attached,wtm,ic1,angles);

  sort_angles(10,angles);
  printf("   angles:");
  for (int i=0;i<10;i++)
    printf(" %5.1f",angles[i]);
  printf("\n");

  if (angles[0]>ANGLE_LINEAR && angles[1]>ANGLE_LINEAR) //square planar + 1
    type = 1;
  else            //not square planar with a hat
    type = 0;

  delete [] angles;
  delete [] attached;

  return type;
}

int ZStruct::find_6c(int wtm, int stm, ICoord ic1)
{
  int type = 0;

  int* attached = new int[12];
  int nattached = find_attached(wtm,stm,ic1,attached);

  double* angles = new double[15];
  find_6c_angles(nattached,attached,wtm,ic1,angles);

  sort_angles(15,angles);
  printf("   angles:");
  for (int i=0;i<15;i++)
    printf(" %5.1f",angles[i]);
  printf("\n");

  if (angles[0]>ANGLE_LINEAR && angles[1]>ANGLE_LINEAR && angles[2]>ANGLE_LINEAR) //octahedral
    type = 1;
  else            //not octahedral
    type = 0;

  delete [] angles;
  delete [] attached;

  return type;
}


void ZStruct::find_6c_angles(int nattached, int* attached, int wtm, ICoord ic1, double* angles)
{
  if (nattached==6)
  {
    int nang = 0;
    for (int i=0;i<6;i++)
    for (int j=0;j<i;j++)
      angles[nang++] = ic1.angle_val(attached[i],wtm,attached[j]);
    printf(" found %i angles \n",nang);
  }
  else if (nattached==7)
  {
    printf("   need to handle single eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2 = attached[2];
    int a3 = attached[3];
    int a4 = attached[4];
    int a5 = attached[5];
    int a6 = attached[5];

    angles[0] = ic1.angle_val(a2,wtm,a3);
    angles[1] = ic1.angle_val(a2,wtm,a4);
    angles[2] = ic1.angle_val(a2,wtm,a5);
    angles[3] = ic1.angle_val(a2,wtm,a6);
    angles[4] = ic1.angle_val(a3,wtm,a4);
    angles[5] = ic1.angle_val(a3,wtm,a5);
    angles[6] = ic1.angle_val(a3,wtm,a6);
    angles[7] = ic1.angle_val(a4,wtm,a5);
    angles[8] = ic1.angle_val(a4,wtm,a6);
    angles[9] = ic1.angle_val(a5,wtm,a6);

    angles[10] = ic1.angle_val_eta2(a2,wtm,a1a,a1b);
    angles[11] = ic1.angle_val_eta2(a3,wtm,a1a,a1b);
    angles[12] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[13] = ic1.angle_val_eta2(a5,wtm,a1a,a1b);
    angles[14] = ic1.angle_val_eta2(a6,wtm,a1a,a1b);
  }
  else if (nattached==8)
  {
    printf("   need to handle double eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];
    int a3 = attached[4];
    int a4 = attached[5];
    int a5 = attached[6];
    int a6 = attached[7];

    angles[0] = ic1.angle_val(a3,wtm,a4);
    angles[1] = ic1.angle_val(a3,wtm,a5);
    angles[2] = ic1.angle_val(a3,wtm,a6);
    angles[3] = ic1.angle_val(a4,wtm,a5);
    angles[4] = ic1.angle_val(a4,wtm,a6);
    angles[5] = ic1.angle_val(a5,wtm,a6);

    angles[6] = ic1.angle_val_eta2(a3,wtm,a1a,a1b);
    angles[7] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[8] = ic1.angle_val_eta2(a5,wtm,a1a,a1b);
    angles[9] = ic1.angle_val_eta2(a6,wtm,a1a,a1b);
    angles[10] = ic1.angle_val_eta2(a3,wtm,a2a,a2b);
    angles[11] = ic1.angle_val_eta2(a4,wtm,a2a,a2b);
    angles[12] = ic1.angle_val_eta2(a5,wtm,a2a,a2b);
    angles[13] = ic1.angle_val_eta2(a6,wtm,a2a,a2b);
    angles[14] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
  }
  else if (nattached==9)
  {
    printf("   need to handle triple eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];
    int a3a = attached[4];
    int a3b = attached[5];
    int a4 = attached[6];
    int a5 = attached[7];
    int a6 = attached[8];

    angles[0] = ic1.angle_val(a4,wtm,a5);
    angles[1] = ic1.angle_val(a4,wtm,a6);
    angles[2] = ic1.angle_val(a5,wtm,a6);

    angles[3] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[4] = ic1.angle_val_eta2(a5,wtm,a1a,a1b);
    angles[5] = ic1.angle_val_eta2(a6,wtm,a1a,a1b);
    angles[6] = ic1.angle_val_eta2(a4,wtm,a2a,a2b);
    angles[7] = ic1.angle_val_eta2(a5,wtm,a2a,a2b);
    angles[8] = ic1.angle_val_eta2(a6,wtm,a2a,a2b);
    angles[9] = ic1.angle_val_eta2(a4,wtm,a3a,a3b);
    angles[10] = ic1.angle_val_eta2(a5,wtm,a3a,a3b);
    angles[11] = ic1.angle_val_eta2(a6,wtm,a3a,a3b);
    angles[12] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
    angles[13] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a3a,a3b);
    angles[14] = ic1.angle_val_eta2_eta2(a2a,a2b,wtm,a3a,a3b);
  }
  else 
  {
    printf("\n ERROR: not yet handling quadruple or higher eta2 \n");
    exit(1);
  }

  return;
}


void ZStruct::find_5c_angles(int nattached, int* attached, int wtm, ICoord ic1, double* angles)
{
  if (nattached==5)
  {
    angles[0] = ic1.angle_val(attached[0],wtm,attached[1]);
    angles[1] = ic1.angle_val(attached[0],wtm,attached[2]);
    angles[2] = ic1.angle_val(attached[0],wtm,attached[3]);
    angles[3] = ic1.angle_val(attached[0],wtm,attached[4]);
    angles[4] = ic1.angle_val(attached[1],wtm,attached[2]);
    angles[5] = ic1.angle_val(attached[1],wtm,attached[3]);
    angles[6] = ic1.angle_val(attached[1],wtm,attached[4]);
    angles[7] = ic1.angle_val(attached[2],wtm,attached[3]);
    angles[8] = ic1.angle_val(attached[2],wtm,attached[4]);
    angles[9] = ic1.angle_val(attached[3],wtm,attached[4]);
  }
  else if (nattached==6)
  {
    printf("   need to handle single eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2 = attached[2];
    int a3 = attached[3];
    int a4 = attached[4];
    int a5 = attached[5];

    angles[0] = ic1.angle_val(a2,wtm,a3);
    angles[1] = ic1.angle_val(a2,wtm,a4);
    angles[2] = ic1.angle_val(a2,wtm,a5);
    angles[3] = ic1.angle_val(a3,wtm,a4);
    angles[4] = ic1.angle_val(a3,wtm,a5);
    angles[5] = ic1.angle_val(a4,wtm,a5);
    angles[6] = ic1.angle_val_eta2(a2,wtm,a1a,a1b);
    angles[7] = ic1.angle_val_eta2(a3,wtm,a1a,a1b);
    angles[8] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[9] = ic1.angle_val_eta2(a5,wtm,a1a,a1b);
  }
  else if (nattached==7)
  {
    printf("   need to handle double eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];
    int a3 = attached[4];
    int a4 = attached[5];
    int a5 = attached[6];

    angles[0] = ic1.angle_val(a3,wtm,a4);
    angles[1] = ic1.angle_val(a3,wtm,a5);
    angles[2] = ic1.angle_val(a4,wtm,a5);
    angles[3] = ic1.angle_val_eta2(a3,wtm,a1a,a1b);
    angles[4] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[5] = ic1.angle_val_eta2(a5,wtm,a1a,a1b);
    angles[6] = ic1.angle_val_eta2(a3,wtm,a2a,a2b);
    angles[7] = ic1.angle_val_eta2(a4,wtm,a2a,a2b);
    angles[8] = ic1.angle_val_eta2(a5,wtm,a2a,a2b);
    angles[9] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
  }
  else if (nattached==8)
  {
    printf("   need to handle triple eta2 \n");
    int a1a = attached[0];
    int a1b = attached[1];
    int a2a = attached[2];
    int a2b = attached[3];
    int a3a = attached[4];
    int a3b = attached[5];
    int a4 = attached[6];
    int a5 = attached[7];

    angles[0] = ic1.angle_val(a4,wtm,a5);
    angles[1] = ic1.angle_val_eta2(a4,wtm,a1a,a1b);
    angles[2] = ic1.angle_val_eta2(a5,wtm,a1a,a1b);
    angles[3] = ic1.angle_val_eta2(a4,wtm,a2a,a2b);
    angles[4] = ic1.angle_val_eta2(a5,wtm,a2a,a2b);
    angles[5] = ic1.angle_val_eta2(a4,wtm,a3a,a3b);
    angles[6] = ic1.angle_val_eta2(a5,wtm,a3a,a3b);
    angles[7] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a2a,a2b);
    angles[8] = ic1.angle_val_eta2_eta2(a1a,a1b,wtm,a3a,a3b);
    angles[9] = ic1.angle_val_eta2_eta2(a2a,a2b,wtm,a3a,a3b);
  }
  else 
  {
    printf("\n ERROR: not yet handling quadruple or higher eta2 \n");
    exit(1);
  }

  return;
}



void ZStruct::sort_angles(int nangles, double* angles)
{
  double* sort = new double[nangles];
  for (int i=0;i<nangles;i++)
    sort[i] = -1000.;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i+1;j++)
    if (sort[j]<angles[i])
    {
      for (int k=nangles-1;k>j;k--)
        sort[k] = sort[k-1];
      sort[j] = angles[i];
      break;
    }
  }
  for (int i=0;i<nangles;i++)
    angles[i] = sort[i];

#if 0
  printf(" sorted angles:");
  for (int i=0;i<nangles;i++)
    printf(" %4.1f",angles[i]);
  printf("\n");
#endif

  delete [] sort;
  return;
}


int ZStruct::find_attached(int wat, int stm, ICoord ic1, int* attached)
{
  int nattached = 0;

  int nbtm = ic1.nbondstm;
  int** btm = ic1.bondstm;
  if (stm==2)
  {
    nbtm = ic1.nbondstm2;
    btm = ic1.bondstm2;
  }

#if 1
 //eta support
  for (int i=0;i<nbtm;i++)
  for (int j=1;j<6;j++)
  if (btm[i][j]>-1)
    attached[nattached++] = btm[i][j];
#else
  for (int i=0;i<ic1.nbonds;i++)
  {
    if (ic1.bonds[i][0]==wat)
      attached[nattached++] = ic1.bonds[i][1];
    else if (ic1.bonds[i][1]==wat)
      attached[nattached++] = ic1.bonds[i][0];
  }
#endif

#if 0
  printf("   attached:");
  for (int i=0;i<nattached;i++)
    printf(" %2i",attached[i]+1);
  printf("\n");
#endif

  return nattached;
}

int ZStruct::determine_geomtype(int wtm, int stm, int coordntm, ICoord ic1)
{
  int geomtype = -1;

  if (coordntm==2)
  {
    geomtype = find_2c(wtm,stm,ic1);
  }
  else if (coordntm==3)
  {
    geomtype = find_3c(wtm,stm,ic1);
  }
  else if (coordntm==4)
  {
    geomtype = find_4c(wtm,stm,ic1);
  }
  else if (coordntm==5)
  {
    geomtype = find_5c(wtm,stm,ic1);
  }
  else if (coordntm==6)
  {
    geomtype = find_6c(wtm,stm,ic1);
  }
  else if (coordntm>6)
  {
    printf("\n ERROR: not yet implemented: coordntm > 6 \n");
    exit(1);
    return 0;
  }

  return geomtype;
}



int ZStruct::get_wtm_addbrktm(int nadd, int* adds, int nbrks, int* brks, ICoord ic1, int& naddtm, int* addtm, int& nbrktm, int* brktm)
{
  int wtm = -1;

  for (int i=0;i<nadd;i++)
  if (ic1.isTM(adds[2*i+0]) || ic1.isTM(adds[2*i+1]))
  {
    int a1 = adds[2*i+0];
    int a2 = adds[2*i+1];
    if (ic1.isTM(a1))
    {
      wtm = a1;
      addtm[2*naddtm+0] = a1;
      addtm[2*naddtm++ +1] = a2;
    }
    else if (ic1.isTM(a2))
    {
      wtm = a2;
      addtm[2*naddtm+0] = a1;
      addtm[2*naddtm++ +1] = a2;
    }
  }
  for (int i=0;i<nbrks;i++)
  if (ic1.isTM(brks[2*i+0]) || ic1.isTM(brks[2*i+1]))
  {
    int b1 = brks[2*i+0];
    int b2 = brks[2*i+1];
    if (ic1.isTM(b1))
    {
      wtm = b1;
      brktm[2*nbrktm+0] = b1;
      brktm[2*nbrktm++ +1] = b2;
    }
    else if (ic1.isTM(b2))
    {
      wtm = b2;
      brktm[2*nbrktm+0] = b1;
      brktm[2*nbrktm++ +1] = b2;
    }
  }

  if (wtm==-1)
  {
    printf("\n failed to find TM, exiting \n");
    exit(1);
  }

  return wtm;
}



#if 0
  else if (nadd==2 && nbrks==1)
  {
    //find remaining group as anchor
    int a1 = ic1.bondstm[0][1];
    if (brks[0]==a1 || brks[1]==a1) //angle must be with remaining group
      a1 = ic1.bondstm[1][1];
    int a2a = adds[0];
    if (adds[0]==wtm)
      a2a = adds[1];
    int a2b = adds[2];
    if (adds[2]==wtm)
      a2b = adds[3];

    ntor = 0;
    nangles = 2;
    angles[0][0] = a1;
    angles[0][1] = wtm;
    angles[0][2] = a2a;
    angles[1][0] = a1;
    angles[1][1] = wtm;
    angles[1][2] = a2b;

    anglev[0] = 120.; //creates Y
    anglev[1] = 120.; 

   //vectoring from wtm pointing to remaining group, a1
    get_vec(a1,wtm,ic1.coords,aprv);

    nfound = 1;
  }
#endif


void ZStruct::find_linear_pairs(int& a1o, int& a2o, int& a3o, int& a4o, int a1, int a2, int a3, int a4, int wtm, ICoord ic1)
{
  double ang1 = ic1.angle_val(a1,wtm,a2);
  double ang2 = ic1.angle_val(a1,wtm,a3);
  double ang3 = ic1.angle_val(a1,wtm,a4);
  double ang4 = ic1.angle_val(a2,wtm,a3);
  double ang5 = ic1.angle_val(a2,wtm,a4);
  double ang6 = ic1.angle_val(a3,wtm,a4);

  if (ang1>ang2 && ang1>ang3)
    a1o = a2;
  else if (ang2>ang1 && ang2>ang3)
    a1o = a3;
  else
    a1o = a4;

  if (ang1>ang4 && ang1>ang5)
    a2o = a1;
  else if (ang4>ang1 && ang4>ang5)
    a2o = a3;
  else
    a2o = a4;

  if (ang2>ang4 && ang2>ang6)
    a3o = a1;
  else if (ang4>ang2 && ang4>ang6)
    a3o = a2;
  else
    a3o = a4;

  if (ang3>ang5 && ang3>ang6)
    a4o = a1;
  else if (ang5>ang3 && ang5>ang6)
    a4o = a2;
  else
    a4o = a3;

  printf("   printing angular opposites: %i-%i %i-%i %i-%i %i-%i \n",a1+1,a1o+1,a2+1,a2o+1,a3+1,a3o+1,a4+1,a4o+1);

  return;
}

void ZStruct::find_oct_planes(int wtm, ICoord ic1, int& p11, int& p12, int& p13, int& p14, int& p21, int& p22, int& p23, int& p24, int& p31, int& p32, int& p33, int& p34)
{
  printf("  in find_oct_planes() \n");
  printf("   NOT IMPLEMENTED \n");
  exit(1);

  p11=p12=p13=p14 = -1;
  p21=p22=p23=p24 = -1;
  p31=p32=p33=p34 = -1;

  return;
}

int ZStruct::find_on_top_5c(int wtm, int stm, ICoord ic1)
{
  int ontop;

  int* attached = new int[10];
  int nattached = find_attached(wtm,stm,ic1,attached);
 
  double* angles = new double[10];
  find_5c_angles(nattached,attached,wtm,ic1,angles);


  if      (angles[0]<ANGLE_LINEAR && angles[1]<ANGLE_LINEAR && angles[2]<ANGLE_LINEAR && angles[3]<ANGLE_LINEAR)
    ontop = ic1.bondstm[0][1];
  else if (angles[0]<ANGLE_LINEAR && angles[4]<ANGLE_LINEAR && angles[5]<ANGLE_LINEAR && angles[6]<ANGLE_LINEAR)
    ontop = ic1.bondstm[1][1];
  else if (angles[1]<ANGLE_LINEAR && angles[4]<ANGLE_LINEAR && angles[7]<ANGLE_LINEAR && angles[8]<ANGLE_LINEAR)
    ontop = ic1.bondstm[2][1];
  else if (angles[2]<ANGLE_LINEAR && angles[5]<ANGLE_LINEAR && angles[7]<ANGLE_LINEAR && angles[9]<ANGLE_LINEAR)
    ontop = ic1.bondstm[3][1];
  else if (angles[3]<ANGLE_LINEAR && angles[6]<ANGLE_LINEAR && angles[8]<ANGLE_LINEAR && angles[9]<ANGLE_LINEAR)
    ontop = ic1.bondstm[4][1];
  else
  {
    printf("  ERROR: couldn't find top atom \n");
    exit(1);
  }


  delete [] angles;
  delete [] attached;

  return ontop;
}

