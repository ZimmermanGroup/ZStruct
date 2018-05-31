// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "align.h"
using namespace std;

//to do
// implement point out for planes


void Align::shuttle_align(int nadd1, int* add1)
{
  if (inited!=3)
  {
    printf(" in shuttle_align, add_third hasn't been called \n");
    exit(1);
  }

#if 0
  printf("\n in shuttle_align, adding:");
  for (int i=0;i<nadd1;i++)
    printf(" %i-%i",add1[2*i+0]+1,add1[2*i+1]+1);
  printf("\n");
  //printf(" in shuttle_align natoms1: %i natoms2: %i \n",natoms1,natoms2);
#endif
#if 0
  printf(" previous geom: \n");
  print_xyz();
#endif

  int na12 = natoms1 + natoms2;

  ///get v1's for all add atoms
  double* v1a = new double[3*nadd1];
  for (int i=0;i<3*nadd1;i++) v1a[i] = 0.;
  double* v1b = new double[3*nadd1];
  for (int i=0;i<3*nadd1;i++) v1b[i] = 0.;
  double* c1 = new double[6]; //center points

  string* anames12 = new string[na12+1];
  int* anumbers12 = new int[na12+1];
  for (int i=0;i<natoms1;i++)
    anames12[i] = anames1[i];
  for (int i=0;i<natoms2;i++)
    anames12[natoms1+i] = anames2[i];
  for (int i=0;i<natoms1;i++)
    anumbers12[i] = anumbers1[i];
  for (int i=0;i<natoms2;i++)
    anumbers12[natoms1+i] = anumbers2[i];

  int* bonded1 = new int[8];
  int* bonded2 = new int[8];
  int* bonded3 = new int[8];

  ICoord ic1,ic2;
  ic1.init(na12,anames12,anumbers12,xyza);
  ic1.make_frags();
  ic1.bond_frags();
  //ic1.ic_create_nobonds();
  ic2.init(natoms3,anames3,anumbers3,xyz3);
  ic2.make_frags();
  ic2.bond_frags();
  //ic2.ic_create_nobonds();

  int npf1 = 0;
  int npf2 = 0;
  int nvf1 = 0;
  int nvf2 = 0;
  int found = 0;
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }

    int same_frag = check_frag_3(a1,a2);
    //printf(" a1: %i a2: %i same_frag: %i \n",a1+1,a2+1,same_frag);
    a2 -= na12;
    if (!same_frag)
    {
      found = 1;
      int nbondsa1, nbondsa2;
      nbondsa1 = get_bonds(a1,ic1,bonded1);
      nbondsa2 = get_bonds(a2,ic2,bonded2);

     //first atom's vector
      for (int j=0;j<nbondsa1;j++)
      {
        v1a[3*nvf1+0] += xyza[3*a1+0] - xyza[3*bonded1[j]+0];
        v1a[3*nvf1+1] += xyza[3*a1+1] - xyza[3*bonded1[j]+1];
        v1a[3*nvf1+2] += xyza[3*a1+2] - xyza[3*bonded1[j]+2];
      }
      if (nbondsa1==2)
      {
         double anglev = ic1.angle_val(bonded1[0],a1,bonded1[1]);
         if (anglev>175.)
         {
           //printf("  using ool v1 \n");
           linear_right(&v1a[3*nvf1],a1,bonded1,xyza);
         }
      }
      if (nbondsa1==3)
      {          
        double imptorv = ic1.torsion_val(bonded1[0],a1,bonded1[1],bonded1[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 \n");
          planar_cross(&v1a[3*nvf1],a1,bonded1,xyza);
          npf1++;
        }
      }
      if (npf1>1)
      {
        align_v1(nvf1+1,v1a);
      }
      double n1 = norm(&v1a[3*nvf1],3);
      for (int j=0;j<3;j++)
        v1a[3*nvf1+j] = v1a[3*nvf1+j] / n1;
#if 0
      printf(" v1a(a1): \n");
      for (int j=0;j<3;j++)
        printf(" %4.3f",v1a[3*nvf1+j]);
      printf("\n");
#endif

      nvf1++;


     //second atom's vector
      for (int j=0;j<nbondsa2;j++)
      {
        v1b[3*nvf2+0] += xyz3[3*a2+0] - xyz3[3*bonded2[j]+0];
        v1b[3*nvf2+1] += xyz3[3*a2+1] - xyz3[3*bonded2[j]+1];
        v1b[3*nvf2+2] += xyz3[3*a2+2] - xyz3[3*bonded2[j]+2];
      }
      if (nbondsa2==2)
      {
         double anglev = ic2.angle_val(bonded2[0],a2,bonded2[1]);
         if (anglev>175.)
         {
           //printf("  using ool v1 (linear angle found) \n");
           linear_right(&v1b[3*nvf2],a2,bonded2,xyz3);
         }
      }
      if (nbondsa2==3)
      {
        double imptorv = ic2.torsion_val(bonded2[0],a2,bonded2[1],bonded2[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 (planar atom found) \n");
          planar_cross(&v1b[3*nvf2],a2,bonded2,xyz3);
          npf2++;
        }
      }
      if (npf2>1)
      {
        align_v1(nvf2+1,v1b);
      }
      n1 = norm(&v1b[3*nvf2],3);
      for (int j=0;j<3;j++)
        v1b[3*nvf2+j] = v1b[3*nvf2+j] / n1;
#if 0
      printf(" v1b(a2): \n");
      for (int j=0;j<3;j++)
        printf(" %4.3f",v1b[3*nvf2+j]);
      printf("\n");
#endif
 
      nvf2++;

    } // if !same_frag
  } //loop i over nadd

  if (!found)
  {
    //printf("  couldn't find anything to align \n");
    return;
  }

#if 0
 //CPMZ NOT IMPLEMENTED
  if (npf1>0)
    point_out(v1a,natoms1,xyz1);
  if (npf2>0)
    point_out(v1b,natoms2,xyz2);
#endif


 //averaging over v1 to get v2
  double* v2 = new double[6];
  for (int i=0;i<6;i++) v2[i] = 0.;
  for (int i=0;i<nadd1;i++)
  {
    for (int j=0;j<3;j++)
      v2[j]   += v1a[3*i+j];
    for (int j=0;j<3;j++)
      v2[3+j] += v1b[3*i+j];
  }

  double n1 = norm(v2,3);
  double n2 = norm(&v2[3],3);
  if (n1<0.000001)
  {
    int a1 = add1[2*(nadd1-1)+0];
    int a2 = add1[2*(nadd1-1)+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    int nbondsa1 = get_bonds(a1,ic1,bonded1);
    linear_right(v2,a1,bonded1,xyza);
    n1 = norm(v2,3);
//    printf(" new norm: %8.6f \n",n1);
  }
  if (n2<0.000001)
  {
    int a1 = add1[2*(nadd1-1)+0];
    int a2 = add1[2*(nadd1-1)+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    a2 -= na12;
    int nbondsa2 = get_bonds(a2,ic2,bonded2);
    linear_right(&v2[3],a2,bonded2,xyz3);
    n2 = norm(&v2[3],3);
//    printf(" new norm: %8.6f \n",n2);
  }
  for (int i=0;i<3;i++)
    v2[i] = v2[i] / n1;
  for (int i=3;i<6;i++)
    v2[i] = v2[i] / n2;

#if 0
  printf("  v2(1):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[j]);
  printf("\n");
  printf("  v2(2):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[3+j]);
  printf("\n");
#endif

  //get central location for each
  int cfound = 0;
  for (int i=0;i<6;i++) c1[i] = 0.;
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    int same_frag = check_frag_3(a1,a2);

    a2 -= na12;
    if (!same_frag)
    {
      //printf(" c1 for %i-%i \n",a1+1,a2+1);
      cfound++;
      for (int j=0;j<3;j++)
        c1[j]   += xyza[3*a1+j];
      for (int j=0;j<3;j++)
        c1[3+j] += xyz3[3*a2+j];
    }
  } //loop i over nadd
  if (cfound>1)
  for (int j=0;j<6;j++)
    c1[j] = c1[j] / cfound;

#if 0
  printf("  c1:");
  for (int j=0;j<6;j++)
    printf(" %4.3f",c1[j]);
  printf("\n");
#endif


  //move geometry to center
  double* xyz12a = new double[3*(na12+2)];
  double* xyz3a = new double[3*(natoms3+2)];
  string* anames12a = new string[na12+2];
  string* anames3a = new string[natoms3+2];
  for (int i=0;i<natoms1;i++) anames12a[i] = anames1[i];
  for (int i=0;i<natoms2;i++) anames12a[natoms1+i] = anames2[i];
  anames12a[na12] = "X";
  anames12a[na12+1] = "X";
  for (int i=0;i<natoms3;i++) anames3a[i] = anames3[i];
  anames3a[natoms3] = "X";
  anames3a[natoms3+1] = "X";

  for (int i=0;i<natoms1;i++)
  {
    for (int j=0;j<3;j++)
      xyz12a[3*i+j] = xyz1[3*i+j] - c1[j];
  }
  for (int i=0;i<natoms2;i++)
  {
    for (int j=0;j<3;j++)
      xyz12a[3*natoms1+3*i+j] = xyz2[3*i+j] - c1[j];
  }
  xyz12a[3*(na12)+0] = 0.; xyz12a[3*(na12)+1] = 0.;
  xyz12a[3*(na12)+2] = 0.; xyz12a[3*(na12+1)+0] = v2[0];
  xyz12a[3*(na12+1)+1] = v2[1]; xyz12a[3*(na12+1)+2] = v2[2];
  //print_xyz_gen(na12+2,anames12a,xyz12a);

  for (int i=0;i<natoms3;i++)
  {
    for (int j=0;j<3;j++)
      xyz3a[3*i+j] = xyz3[3*i+j] - c1[3+j];
  }
  xyz3a[3*(natoms3)+0] = 0.; xyz3a[3*(natoms3)+1] = 0.;
  xyz3a[3*(natoms3)+2] = 0.; xyz3a[3*(natoms3+1)+0] = v2[3];
  xyz3a[3*(natoms3+1)+1] = v2[4]; xyz3a[3*(natoms3+1)+2] = v2[5];
  //print_xyz_gen(natoms3+2,anames3a,xyz3a);

  //align v2 along x/-x direction
  int t1 = na12;
  int t2 = na12+1;
  align_to_x(na12+2,t1,t2,xyz12a,anames12a,1,10.);
  t1 = natoms3;
  t2 = natoms3+1;
  align_to_x(natoms3+2,t1,t2,xyz3a,anames3a,-1,10.);

  //print_xyz_gen(na12+2,anames12a,xyz12a);
  //print_xyz_gen(natoms3+2,anames3a,xyz3a);

  //face each other at X Angstroms
  double X = 8.;
  for (int i=0;i<3*natoms1;i++)
    xyz1[i] = xyz12a[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2[i] = xyz12a[3*natoms1+i];
  for (int i=0;i<3*natoms3;i++)
    xyz3[i] = xyz3a[i];
  for (int i=0;i<natoms3;i++)
    xyz3[3*i+0] += X;


 //prepare for vdw opt
  string* anamesp = new string[na12+natoms3];
  int* anumbersp = new int[na12+natoms3];
  for (int i=0;i<natoms1;i++)
    anamesp[i] = anames1[i];
  for (int i=0;i<natoms2;i++)
    anamesp[natoms1+i] = anames2[i];
  for (int i=0;i<natoms3;i++)
    anamesp[na12+i] = anames3[i];
  for (int i=0;i<natoms1;i++)
    anumbersp[i] = anumbers1[i];
  for (int i=0;i<natoms2;i++)
    anumbersp[natoms1+i] = anumbers2[i];
  for (int i=0;i<natoms3;i++)
    anumbersp[na12+i] = anumbers3[i];
 //full pair geometry, before vdw opt
  for (int i=0;i<3*natoms1;i++)
    xyza3[i] = xyz1[i];
  for (int i=0;i<3*natoms2;i++)
    xyza3[3*natoms1+i] = xyz2[i];
  for (int i=0;i<3*natoms3;i++)
    xyza3[3*(na12)+i] = xyz3[i];


  ICoord icp;
  icp.init(na12+natoms3,anamesp,anumbersp,xyza3);

  //print_xyz();

  //rotate about central axis
  if (nadd1>1)
  {
    int a1 = add1[0];
    int a2 = add1[1];
    int tmp1 = a1; if (a1>a2) { a1 = a2; a2 = tmp1; }
    int a3 = add1[2];
    int a4 = add1[3];
    int tmp2 = a3; if (a3>a4) { a3 = a4; a4 = tmp2; }

    int same_frag = check_frag_3(a1,a2) + check_frag_3(a3,a4);
    if (!same_frag && a1!=a3 && a2!=a4)
    {
      double torv = icp.torsion_val(a1,a3,a4,a2);
      //printf("  central torv: %4.2f (%i %i %i %i) \n",torv,a1+1,a3+1,a4+1,a2+1);
      rotate_around_x(na12,natoms3,-torv,icp.coords);
    }
  }

  double* vx = new double[3];
  vx[1] = vx[2] = 0.;
  vx[0] = 1.0;
  vdw_vector_opt(na12,natoms3,vx,icp);
  delete [] vx;

  for (int i=0;i<3*(na12+natoms3);i++)
    xyza3[i] = icp.coords[i];

  for (int i=0;i<3*natoms1;i++)
    xyz1[i] = xyza3[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2[i] = xyza3[3*natoms1+i];
  for (int i=0;i<3*natoms3;i++)
    xyz3[i] = xyza3[3*(na12)+i];

  //print_xyz_gen(na12,icp.anames,xyza);

  int badgeom = 0;
  for (int i=0;i<3*natoms1;i++)
  if (xyz1[i]!=xyz1[i])
  {
    badgeom = 1;
    break;
  }
  for (int i=0;i<3*natoms2;i++)
  if (xyz2[i]!=xyz2[i])
  {
    badgeom = 1;
    break;
  }
  for (int i=0;i<3*natoms3;i++)
  if (xyz3[i]!=xyz3[i])
  {
    badgeom = 1;
    break;
  }

  if (badgeom)
  {
    printf(" ERROR: nan detected in shuttle_align, exiting! \n");
    print_xyz_gen(natoms1,anames1,xyz1);
    print_xyz_gen(natoms2,anames2,xyz2);
    print_xyz_gen(natoms3,anames3,xyz3);
    exit(1);
  }

  delete [] anames12;
  delete [] anumbers12;

  delete [] anames12a;
  delete [] anames3a;
  delete [] xyz12a;
  delete [] xyz3a;

  delete [] v2;
  delete [] v1a;
  delete [] v1b;
  delete [] bonded1;
  delete [] bonded2;
  delete [] bonded3;

  ic1.freemem();
  ic2.freemem();

  delete [] anamesp;
  delete [] anumbersp;
  icp.freemem();


  return;
}


void Align::add_third(int natoms3i, string* anames3i, int* anumbers3i, double* xyz3i)
{
  if (!inited)
  {
    printf(" WARNING: cannot add_third without initing first \n");
    return;
  }
  if (inited==3)
  {
    if (anames3!=NULL) delete [] anames3;
    if (anumbers3!=NULL) delete [] anumbers3;
    if (xyz3!=NULL) delete [] xyz3;
    if (xyza3!=NULL) delete [] xyza3;
  }
  inited = 3;

  natoms3 = natoms3i;
  anames3 = new string[natoms3];
  anumbers3 = new int[natoms3];
  xyz3 = new double[3*natoms3];
  xyza3 = new double[3*(natoms1+natoms2+natoms3)];

  for (int i=0;i<natoms3;i++)
  {
    anames3[i] = anames3i[i];
    anumbers3[i] = anumbers3i[i];
  }
  for (int i=0;i<3*natoms3;i++)
    xyz3[i] = xyz3i[i];

  return;
}

void Align::vdw_vector_opt(int natoms1, int natoms2, double* v1, ICoord icp)
{
  //printf(" in vdw_vector_opt \n");

  double E = 0.;
  double pE = 0.;
  double THRESH = 0.001;
  double emin = 1000000.;
  for (int i=0;i<200;i++)
  {
    pE = E;
    E = icp.mm_energy();
    //printf(" E: %6.5f \n",E);

    if (E<emin) emin = E;
    if (E>pE+THRESH && E>emin) break;

    for (int j=0;j<natoms2;j++)
    for (int k=0;k<3;k++)
      icp.coords[3*natoms1+3*j+k] -= 0.1*v1[k];

   // print_xyz_gen(icp.natoms,icp.anames,icp.coords);
  }

  return;
}

void Align::point_out(double* v1, int natoms, double* xyz)
{

  return;
}


void Align::align_zero()
{
  if (inited==0) return;

  //printf(" in align_zero \n");

 //get vector pointed along COM's
 //move close


 //full pair geometry
  for (int i=0;i<3*natoms1;i++)
    xyza[i] = xyz1[i];
  for (int i=0;i<3*natoms2;i++)
    xyza[3*natoms1+i] = xyz2[i];

  return;
}



int Align::add_align_v(int nadd1, int* add1, int wtm, double* aprv)
{
  if (inited==0)
  {
    printf(" WARNING: in add_align, not init'd \n");
    return 0;
  }
 //revert to initial geometries
  for (int i=0;i<3*natoms1;i++)
    xyz1[i] = xyz1s[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2[i] = xyz2s[i];

#if 1
  printf("   in add_align_v, adding:");
  for (int i=0;i<nadd1;i++)
    printf(" %i-%i",add1[2*i+0]+1,add1[2*i+1]+1);
  printf("\n");
  //printf(" in add_align natoms1: %i natoms2: %i \n",natoms1,natoms2);
#endif

  ///get v1's for all add atoms
  double* v1a = new double[3*nadd1];
  for (int i=0;i<3*nadd1;i++) v1a[i] = 0.;
  double* v1b = new double[3*nadd1];
  for (int i=0;i<3*nadd1;i++) v1b[i] = 0.;
  double* c1 = new double[6]; //center points

  int* bonded1 = new int[8];
  int* bonded2 = new int[8];
  int* bonded3 = new int[8];

  ICoord ic1,ic2;
  ic1.init(natoms1,anames1,anumbers1,xyz1);
  ic1.make_frags();
  ic1.bond_frags();
  //ic1.ic_create_nobonds();
  //ic1.print_bonds();
  ic2.init(natoms2,anames2,anumbers2,xyz2);
  ic2.make_frags();
  ic2.bond_frags();
  //ic2.ic_create_nobonds();
  //ic2.print_bonds();

  int npf1 = 0;
  int npf2 = 0;
  int nvf1 = 0;
  int nvf2 = 0;
  int found = 0;
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }

    int same_frag = check_frag(a1,a2);
    //printf(" a1: %i a2: %i same_frag: %i \n",a1+1,a2+1,same_frag);
    a2 -= natoms1;
    if (!same_frag)
    {
      found = 1;
      int nbondsa1, nbondsa2;
      nbondsa1 = get_bonds(a1,ic1,bonded1);
      nbondsa2 = get_bonds(a2,ic2,bonded2);

     //first atom's vector
      for (int j=0;j<nbondsa1;j++)
      {
        v1a[3*nvf1+0] += xyz1[3*a1+0] - xyz1[3*bonded1[j]+0];
        v1a[3*nvf1+1] += xyz1[3*a1+1] - xyz1[3*bonded1[j]+1];
        v1a[3*nvf1+2] += xyz1[3*a1+2] - xyz1[3*bonded1[j]+2];
      }
#if 0
      if (nbondsa1==1 && ic1.coordn[bonded1[0]]==3)
      {
        int a1c = bonded1[0];
        int nbondsa1c = get_bonds(a1c,ic1,bonded3)
        double imptorv = ic1.torsion_val(bonded3[0],a1,bonded3[1],bonded3[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 \n");
          planar_cross(&v1a[3*nvf1],a3,bonded3,xyz1);
        }
      }
#endif
      if (nbondsa1==2)
      {
         double anglev = ic1.angle_val(bonded1[0],a1,bonded1[1]);
         if (anglev>175.)
         {
           //printf("  using ool v1 \n");
           linear_right(&v1a[3*nvf1],a1,bonded1,xyz1);
         }
      }
      if (nbondsa1==3)
      { 
        double imptorv = ic1.torsion_val(bonded1[0],a1,bonded1[1],bonded1[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 \n");
          planar_cross(&v1a[3*nvf1],a1,bonded1,xyz1);
          npf1++;
        }
      }
      if (npf1>1)
      {
        align_v1(nvf1+1,v1a);
      }
      double n1 = norm(&v1a[3*nvf1],3);
      for (int j=0;j<3;j++)
        v1a[3*nvf1+j] = v1a[3*nvf1+j] / n1;
#if 0
      printf(" v1a(a1): \n");
      for (int j=0;j<3;j++)
        printf(" %4.3f",v1a[3*nvf1+j]);
      printf("\n");
#endif

      nvf1++;


     //second atom's vector
      for (int j=0;j<nbondsa2;j++)
      {
        v1b[3*nvf2+0] += xyz2[3*a2+0] - xyz2[3*bonded2[j]+0];
        v1b[3*nvf2+1] += xyz2[3*a2+1] - xyz2[3*bonded2[j]+1];
        v1b[3*nvf2+2] += xyz2[3*a2+2] - xyz2[3*bonded2[j]+2];
      }
      if (nbondsa2==2)
      {
         double anglev = ic2.angle_val(bonded2[0],a2,bonded2[1]);
         if (anglev>175.)
         {
           //printf("  using ool v1 (linear angle found) \n");
           linear_right(&v1b[3*nvf2],a2,bonded2,xyz2);
         }
      }
      if (nbondsa2==3)
      {
        double imptorv = ic2.torsion_val(bonded2[0],a2,bonded2[1],bonded2[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 (planar atom found) \n");
          planar_cross(&v1b[3*nvf2],a2,bonded2,xyz2);
          npf2++;
        }
      }
      if (npf2>1)
      {
        align_v1(nvf2+1,v1b);
      }
      n1 = norm(&v1b[3*nvf2],3);
      for (int j=0;j<3;j++)
        v1b[3*nvf2+j] = v1b[3*nvf2+j] / n1;
#if 0
      printf(" v1b(a2): \n");
      for (int j=0;j<3;j++)
        printf(" %4.3f",v1b[3*nvf2+j]);
      printf("\n");
#endif
 
      nvf2++;

    } // if !same_frag
  } //loop i over nadd

  if (!found)
  {
    //printf("  couldn't find anything to align \n");
    return 0;
  }

#if 0
 //CPMZ NOT IMPLEMENTED
  if (npf1>0)
    point_out(v1a,natoms1,xyz1);
  if (npf2>0)
    point_out(v1b,natoms2,xyz2);
#endif


 //averaging over v1 to get v2
  double* v2 = new double[6];
  for (int i=0;i<6;i++) v2[i] = 0.;
  for (int i=0;i<nadd1;i++)
  {
    for (int j=0;j<3;j++)
      v2[j]   += v1a[3*i+j];
    for (int j=0;j<3;j++)
      v2[3+j] += v1b[3*i+j];
  }

  double n1 = norm(v2,3);
  double n2 = norm(&v2[3],3);
  if (n1<0.000001)
  {
    int a1 = add1[2*(nadd1-1)+0];
    int a2 = add1[2*(nadd1-1)+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    int nbondsa1 = get_bonds(a1,ic1,bonded1);
    linear_right(v2,a1,bonded1,xyz1);
    n1 = norm(v2,3);
//    printf(" new norm: %8.6f \n",n1);
  }
  if (n2<0.000001)
  {
    int a1 = add1[2*(nadd1-1)+0];
    int a2 = add1[2*(nadd1-1)+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    a2 -= natoms1;
    int nbondsa2 = get_bonds(a2,ic2,bonded2);
    linear_right(&v2[3],a2,bonded2,xyz2);
    n2 = norm(&v2[3],3);
//    printf(" new norm: %8.6f \n",n2);
  }
  for (int i=0;i<3;i++)
    v2[i] = v2[i] / n1;
  for (int i=3;i<6;i++)
    v2[i] = v2[i] / n2;

#if 0
  printf("     v2(1):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[j]);
  printf("\n");
  printf("     v2(2):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[3+j]);
  printf("\n");
#endif

  //get central location for each
  int cfound = 0;
  for (int i=0;i<6;i++) c1[i] = 0.;
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    int same_frag = check_frag(a1,a2);

    a2 -= natoms1;
    if (!same_frag)
    {
      //printf(" c1 for %i-%i \n",a1+1,a2+1);
      cfound++;
      for (int j=0;j<3;j++)
        c1[j]   += xyz1[3*a1+j];
      for (int j=0;j<3;j++)
        c1[3+j] += xyz2[3*a2+j];
    }
  } //loop i over nadd
  if (cfound>1)
  for (int j=0;j<6;j++)
    c1[j] = c1[j] / cfound;

#if 0
  printf("  c1:");
  for (int j=0;j<6;j++)
    printf(" %4.3f",c1[j]);
  printf("\n");
#endif


  //find TM fragment, replace c1 and v terms with wtm center and aprv
  int f1 = 0;
  for (int i=0;i<natoms1;i++)
  if (ic1.isTM(i))
    f1 = 1;
  int f2 = 0;
  for (int i=0;i<natoms2;i++)
  if (ic2.isTM(i))
    f2 = 1;

  int TMaddo = 0; //is TM adding to other reactant?
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }

    int same_frag = check_frag(a1,a2);
    //printf(" a1: %i a2: %i same_frag: %i \n",a1+1,a2+1,same_frag);
    a2 -= natoms1;
    if (!same_frag)
    {
      if (ic1.isTM(a1))
      {
        TMaddo++;
        f1++;
      }
      if (ic2.isTM(a2))
      {
        TMaddo++;
        f2++;
      }
    }
  }

  //printf("    fragment with TM: %i %i add across TM-other: %i \n",f1,f2,TMaddo);
  if (f1==f2 && TMaddo)
  {
    printf("   ERROR: neither or both fragments have TMs! \n");
    exit(1);
  }
  else if (f1>0 && TMaddo)
  {
    c1[0] = xyz1[3*wtm+0];
    c1[1] = xyz1[3*wtm+1];
    c1[2] = xyz1[3*wtm+2];
    v2[0] = aprv[0];
    v2[1] = aprv[1];
    v2[2] = aprv[2];
  }
  else if (f2>0 && TMaddo)
  {
    wtm  -= natoms1;
    c1[3] = xyz2[3*wtm+0];
    c1[4] = xyz2[3*wtm+1];
    c1[5] = xyz2[3*wtm+2];
    v2[3] = aprv[0];
    v2[4] = aprv[1];
    v2[5] = aprv[2];
  }
  else
    printf("   not using TM-specific vector \n");

#if 0
  printf("  TM c1:");
  for (int j=0;j<6;j++)
    printf(" %4.3f",c1[j]);
  printf("\n");
#endif
#if 0
  printf("  TM v2(1):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[j]);
  printf("\n");
  printf("  TM v2(2):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[3+j]);
  printf("\n");
#endif

  //move geometry to center
  double* xyz1a = new double[3*(natoms1+2)];
  double* xyz2a = new double[3*(natoms2+2)];
  string* anames1a = new string[natoms1+2];
  string* anames2a = new string[natoms2+2];
  for (int i=0;i<natoms1;i++) anames1a[i] = anames1[i];
  anames1a[natoms1] = "X";
  anames1a[natoms1+1] = "X";
  for (int i=0;i<natoms2;i++) anames2a[i] = anames2[i];
  anames2a[natoms2] = "X";
  anames2a[natoms2+1] = "X";

  for (int i=0;i<natoms1;i++)
  {
    for (int j=0;j<3;j++)
      xyz1a[3*i+j] = xyz1[3*i+j] - c1[j];
  }
  xyz1a[3*(natoms1)+0] = 0.; xyz1a[3*(natoms1)+1] = 0.;
  xyz1a[3*(natoms1)+2] = 0.; xyz1a[3*(natoms1+1)+0] = v2[0];
  xyz1a[3*(natoms1+1)+1] = v2[1]; xyz1a[3*(natoms1+1)+2] = v2[2];
  //print_xyz_gen(natoms1+2,anames1a,xyz1a);

  for (int i=0;i<natoms2;i++)
  {
    for (int j=0;j<3;j++)
      xyz2a[3*i+j] = xyz2[3*i+j] - c1[3+j];
  }
  xyz2a[3*(natoms2)+0] = 0.; xyz2a[3*(natoms2)+1] = 0.;
  xyz2a[3*(natoms2)+2] = 0.; xyz2a[3*(natoms2+1)+0] = v2[3];
  xyz2a[3*(natoms2+1)+1] = v2[4]; xyz2a[3*(natoms2+1)+2] = v2[5];
  //print_xyz_gen(natoms2+2,anames2a,xyz2a);


  //align v2 along x/-x direction
  int t1 = natoms1;
  int t2 = natoms1+1;
  align_to_x(natoms1+2,t1,t2,xyz1a,anames1a,1,10.);
  t1 = natoms2;
  t2 = natoms2+1;
  align_to_x(natoms2+2,t1,t2,xyz2a,anames2a,-1,10.);

  //print_xyz_gen(natoms1+2,anames1a,xyz1a);
  //print_xyz_gen(natoms2+2,anames2a,xyz2a);

  //face each other at X Angstroms
  double X = 18.;
  for (int i=0;i<3*natoms1;i++)
    xyz1[i] = xyz1a[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2[i] = xyz2a[i];
  for (int i=0;i<natoms2;i++)
    xyz2[3*i+0] += X;


 //prepare for vdw opt
  string* anamesp = new string[natoms1+natoms2];
  int* anumbersp = new int[natoms1+natoms2];
  for (int i=0;i<natoms1;i++)
    anamesp[i] = anames1[i];
  for (int i=0;i<natoms2;i++)
    anamesp[natoms1+i] = anames2[i];
  for (int i=0;i<natoms1;i++)
    anumbersp[i] = anumbers1[i];
  for (int i=0;i<natoms2;i++)
    anumbersp[natoms1+i] = anumbers2[i];
 //full pair geometry, before vdw opt
  for (int i=0;i<3*natoms1;i++)
    xyza[i] = xyz1[i];
  for (int i=0;i<3*natoms2;i++)
    xyza[3*natoms1+i] = xyz2[i];

  ICoord icp;
  icp.init(natoms1+natoms2,anamesp,anumbersp,xyza);

  //print_xyz();

  //rotate about central axis
  if (nadd1>1)
  {
    int a1 = add1[0];
    int a2 = add1[1];
    int tmp1 = a1; if (a1>a2) { a1 = a2; a2 = tmp1; }
    int a3 = add1[2];
    int a4 = add1[3];
    int tmp2 = a3; if (a3>a4) { a3 = a4; a4 = tmp2; }

    int same_frag = check_frag(a1,a2) + check_frag(a3,a4);
    if (!same_frag && a1!=a3 && a2!=a4)
    {
      double torv = icp.torsion_val(a1,a3,a4,a2);
      //printf("  central torv: %4.2f (%i %i %i %i) \n",torv,a1+1,a3+1,a4+1,a2+1);
      rotate_around_x(natoms1,natoms2,-torv,icp.coords);
    }
  }

  double* vx = new double[3];
  vx[1] = vx[2] = 0.;
  vx[0] = 1.0;
  vdw_vector_opt(natoms1,natoms2,vx,icp);
  delete [] vx;

  for (int i=0;i<3*(natoms1+natoms2);i++)
    xyza[i] = icp.coords[i];

  for (int i=0;i<3*natoms1;i++)
    xyz1[i] = xyza[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2[i] = xyza[3*natoms1+i];

  //print_xyz_gen(natoms1+natoms2,icp.anames,xyza);

  int badgeom = 0;
  for (int i=0;i<3*natoms1;i++)
  if (xyz1[i]!=xyz1[i])
  {
    badgeom = 1;
    break;
  }
  for (int i=0;i<3*natoms2;i++)
  if (xyz2[i]!=xyz2[i])
  {
    badgeom = 1;
    break;
  }

  if (badgeom)
  {
    printf(" ERROR: nan detected in add_align, exiting! \n");
    print_xyz_gen(natoms1,anames1,xyz1);
    print_xyz_gen(natoms2,anames2,xyz2);
    exit(1);
  }

  delete [] anames1a;
  delete [] anames2a;
  delete [] xyz1a;
  delete [] xyz2a;

  delete [] v2;
  delete [] v1a;
  delete [] v1b;
  delete [] bonded1;
  delete [] bonded2;
  delete [] bonded3;

  ic1.freemem();
  ic2.freemem();

  delete [] anamesp;
  delete [] anumbersp;
  icp.freemem();


  return TMaddo;
}




void Align::add_align(int nadd1, int* add1)
{
  if (inited==0)
  {
    printf(" WARNING: in add_align, not init'd \n");
    return;
  }

#if 0
  printf("\n in add_align, adding:");
  for (int i=0;i<nadd1;i++)
    printf(" %i-%i",add1[2*i+0]+1,add1[2*i+1]+1);
  printf("\n");
  //printf(" in add_align natoms1: %i natoms2: %i \n",natoms1,natoms2);
#endif

  ///get v1's for all add atoms
  double* v1a = new double[3*nadd1];
  for (int i=0;i<3*nadd1;i++) v1a[i] = 0.;
  double* v1b = new double[3*nadd1];
  for (int i=0;i<3*nadd1;i++) v1b[i] = 0.;
  double* c1 = new double[6]; //center points

  int* bonded1 = new int[8];
  int* bonded2 = new int[8];
  int* bonded3 = new int[8];

  ICoord ic1,ic2;
  ic1.init(natoms1,anames1,anumbers1,xyz1);
  ic1.make_frags();
  ic1.bond_frags();
  //ic1.ic_create_nobonds();
  //ic1.print_bonds();
  ic2.init(natoms2,anames2,anumbers2,xyz2);
  ic2.make_frags();
  ic2.bond_frags();
  //ic2.ic_create_nobonds();
  //ic2.print_bonds();

  int npf1 = 0;
  int npf2 = 0;
  int nvf1 = 0;
  int nvf2 = 0;
  int found = 0;
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }

    int same_frag = check_frag(a1,a2);
    //printf(" a1: %i a2: %i same_frag: %i \n",a1+1,a2+1,same_frag);
    a2 -= natoms1;
    if (!same_frag)
    {
      found = 1;
      int nbondsa1, nbondsa2;
      nbondsa1 = get_bonds(a1,ic1,bonded1);
      nbondsa2 = get_bonds(a2,ic2,bonded2);

     //first atom's vector
      for (int j=0;j<nbondsa1;j++)
      {
        v1a[3*nvf1+0] += xyz1[3*a1+0] - xyz1[3*bonded1[j]+0];
        v1a[3*nvf1+1] += xyz1[3*a1+1] - xyz1[3*bonded1[j]+1];
        v1a[3*nvf1+2] += xyz1[3*a1+2] - xyz1[3*bonded1[j]+2];
      }
#if 0
      if (nbondsa1==1 && ic1.coordn[bonded1[0]]==3)
      {
        int a1c = bonded1[0];
        int nbondsa1c = get_bonds(a1c,ic1,bonded3)
        double imptorv = ic1.torsion_val(bonded3[0],a1,bonded3[1],bonded3[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 \n");
          planar_cross(&v1a[3*nvf1],a3,bonded3,xyz1);
        }
      }
#endif
      if (nbondsa1==2)
      {
         double anglev = ic1.angle_val(bonded1[0],a1,bonded1[1]);
         if (anglev>175.)
         {
           //printf("  using ool v1 \n");
           linear_right(&v1a[3*nvf1],a1,bonded1,xyz1);
         }
      }
      if (nbondsa1==3)
      {          
        double imptorv = ic1.torsion_val(bonded1[0],a1,bonded1[1],bonded1[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 \n");
          planar_cross(&v1a[3*nvf1],a1,bonded1,xyz1);
          npf1++;
        }
      }
      if (npf1>1)
      {
        align_v1(nvf1+1,v1a);
      }
      double n1 = norm(&v1a[3*nvf1],3);
      for (int j=0;j<3;j++)
        v1a[3*nvf1+j] = v1a[3*nvf1+j] / n1;
#if 0
      printf(" v1a(a1): \n");
      for (int j=0;j<3;j++)
        printf(" %4.3f",v1a[3*nvf1+j]);
      printf("\n");
#endif

      nvf1++;


     //second atom's vector
      for (int j=0;j<nbondsa2;j++)
      {
        v1b[3*nvf2+0] += xyz2[3*a2+0] - xyz2[3*bonded2[j]+0];
        v1b[3*nvf2+1] += xyz2[3*a2+1] - xyz2[3*bonded2[j]+1];
        v1b[3*nvf2+2] += xyz2[3*a2+2] - xyz2[3*bonded2[j]+2];
      }
      if (nbondsa2==2)
      {
         double anglev = ic2.angle_val(bonded2[0],a2,bonded2[1]);
         if (anglev>175.)
         {
           //printf("  using ool v1 (linear angle found) \n");
           linear_right(&v1b[3*nvf2],a2,bonded2,xyz2);
         }
      }
      if (nbondsa2==3)
      {
        double imptorv = ic2.torsion_val(bonded2[0],a2,bonded2[1],bonded2[2]);
        //printf("  imptorv: %4.2f \n",imptorv);
        if (fabs(imptorv)>175.)
        {
          //printf("  using oop v1 (planar atom found) \n");
          planar_cross(&v1b[3*nvf2],a2,bonded2,xyz2);
          npf2++;
        }
      }
      if (npf2>1)
      {
        align_v1(nvf2+1,v1b);
      }
      n1 = norm(&v1b[3*nvf2],3);
      for (int j=0;j<3;j++)
        v1b[3*nvf2+j] = v1b[3*nvf2+j] / n1;
#if 0
      printf(" v1b(a2): \n");
      for (int j=0;j<3;j++)
        printf(" %4.3f",v1b[3*nvf2+j]);
      printf("\n");
#endif
 
      nvf2++;

    } // if !same_frag
  } //loop i over nadd

  if (!found)
  {
    //printf("  couldn't find anything to align \n");
    return;
  }

#if 0
 //CPMZ NOT IMPLEMENTED
  if (npf1>0)
    point_out(v1a,natoms1,xyz1);
  if (npf2>0)
    point_out(v1b,natoms2,xyz2);
#endif


 //averaging over v1 to get v2
  double* v2 = new double[6];
  for (int i=0;i<6;i++) v2[i] = 0.;
  for (int i=0;i<nadd1;i++)
  {
    for (int j=0;j<3;j++)
      v2[j]   += v1a[3*i+j];
    for (int j=0;j<3;j++)
      v2[3+j] += v1b[3*i+j];
  }

  double n1 = norm(v2,3);
  double n2 = norm(&v2[3],3);
  if (n1<0.000001)
  {
    int a1 = add1[2*(nadd1-1)+0];
    int a2 = add1[2*(nadd1-1)+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    int nbondsa1 = get_bonds(a1,ic1,bonded1);
    linear_right(v2,a1,bonded1,xyz1);
    n1 = norm(v2,3);
//    printf(" new norm: %8.6f \n",n1);
  }
  if (n2<0.000001)
  {
    int a1 = add1[2*(nadd1-1)+0];
    int a2 = add1[2*(nadd1-1)+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    a2 -= natoms1;
    int nbondsa2 = get_bonds(a2,ic2,bonded2);
    linear_right(&v2[3],a2,bonded2,xyz2);
    n2 = norm(&v2[3],3);
//    printf(" new norm: %8.6f \n",n2);
  }
  for (int i=0;i<3;i++)
    v2[i] = v2[i] / n1;
  for (int i=3;i<6;i++)
    v2[i] = v2[i] / n2;

#if 0
  printf("  v2(1):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[j]);
  printf("\n");
  printf("  v2(2):");
  for (int j=0;j<3;j++)
    printf(" %4.3f",v2[3+j]);
  printf("\n");
#endif

  //get central location for each
  int cfound = 0;
  for (int i=0;i<6;i++) c1[i] = 0.;
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int tmp = a1; if (a1>a2) { a1 = a2; a2 = tmp; }
    int same_frag = check_frag(a1,a2);

    a2 -= natoms1;
    if (!same_frag)
    {
      //printf(" c1 for %i-%i \n",a1+1,a2+1);
      cfound++;
      for (int j=0;j<3;j++)
        c1[j]   += xyz1[3*a1+j];
      for (int j=0;j<3;j++)
        c1[3+j] += xyz2[3*a2+j];
    }
  } //loop i over nadd
  if (cfound>1)
  for (int j=0;j<6;j++)
    c1[j] = c1[j] / cfound;

#if 0
  printf("  c1:");
  for (int j=0;j<6;j++)
    printf(" %4.3f",c1[j]);
  printf("\n");
#endif


  //move geometry to center
  double* xyz1a = new double[3*(natoms1+2)];
  double* xyz2a = new double[3*(natoms2+2)];
  string* anames1a = new string[natoms1+2];
  string* anames2a = new string[natoms2+2];
  for (int i=0;i<natoms1;i++) anames1a[i] = anames1[i];
  anames1a[natoms1] = "X";
  anames1a[natoms1+1] = "X";
  for (int i=0;i<natoms2;i++) anames2a[i] = anames2[i];
  anames2a[natoms2] = "X";
  anames2a[natoms2+1] = "X";

  for (int i=0;i<natoms1;i++)
  {
    for (int j=0;j<3;j++)
      xyz1a[3*i+j] = xyz1[3*i+j] - c1[j];
  }
  xyz1a[3*(natoms1)+0] = 0.; xyz1a[3*(natoms1)+1] = 0.;
  xyz1a[3*(natoms1)+2] = 0.; xyz1a[3*(natoms1+1)+0] = v2[0];
  xyz1a[3*(natoms1+1)+1] = v2[1]; xyz1a[3*(natoms1+1)+2] = v2[2];
  //print_xyz_gen(natoms1+2,anames1a,xyz1a);

  for (int i=0;i<natoms2;i++)
  {
    for (int j=0;j<3;j++)
      xyz2a[3*i+j] = xyz2[3*i+j] - c1[3+j];
  }
  xyz2a[3*(natoms2)+0] = 0.; xyz2a[3*(natoms2)+1] = 0.;
  xyz2a[3*(natoms2)+2] = 0.; xyz2a[3*(natoms2+1)+0] = v2[3];
  xyz2a[3*(natoms2+1)+1] = v2[4]; xyz2a[3*(natoms2+1)+2] = v2[5];
  //print_xyz_gen(natoms2+2,anames2a,xyz2a);


  //align v2 along x/-x direction
  int t1 = natoms1;
  int t2 = natoms1+1;
  align_to_x(natoms1+2,t1,t2,xyz1a,anames1a,1,10.);
  t1 = natoms2;
  t2 = natoms2+1;
  align_to_x(natoms2+2,t1,t2,xyz2a,anames2a,-1,10.);

  //print_xyz_gen(natoms1+2,anames1a,xyz1a);
  //print_xyz_gen(natoms2+2,anames2a,xyz2a);

  //face each other at X Angstroms
  double X = 8.;
  for (int i=0;i<3*natoms1;i++)
    xyz1[i] = xyz1a[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2[i] = xyz2a[i];
  for (int i=0;i<natoms2;i++)
    xyz2[3*i+0] += X;


 //prepare for vdw opt
  string* anamesp = new string[natoms1+natoms2];
  int* anumbersp = new int[natoms1+natoms2];
  for (int i=0;i<natoms1;i++)
    anamesp[i] = anames1[i];
  for (int i=0;i<natoms2;i++)
    anamesp[natoms1+i] = anames2[i];
  for (int i=0;i<natoms1;i++)
    anumbersp[i] = anumbers1[i];
  for (int i=0;i<natoms2;i++)
    anumbersp[natoms1+i] = anumbers2[i];
 //full pair geometry, before vdw opt
  for (int i=0;i<3*natoms1;i++)
    xyza[i] = xyz1[i];
  for (int i=0;i<3*natoms2;i++)
    xyza[3*natoms1+i] = xyz2[i];

  ICoord icp;
  icp.init(natoms1+natoms2,anamesp,anumbersp,xyza);

  //print_xyz();

  //rotate about central axis
  if (nadd1>1)
  {
    int a1 = add1[0];
    int a2 = add1[1];
    int tmp1 = a1; if (a1>a2) { a1 = a2; a2 = tmp1; }
    int a3 = add1[2];
    int a4 = add1[3];
    int tmp2 = a3; if (a3>a4) { a3 = a4; a4 = tmp2; }

    int same_frag = check_frag(a1,a2) + check_frag(a3,a4);
    if (!same_frag && a1!=a3 && a2!=a4)
    {
      double torv = icp.torsion_val(a1,a3,a4,a2);
      //printf("  central torv: %4.2f (%i %i %i %i) \n",torv,a1+1,a3+1,a4+1,a2+1);
      rotate_around_x(natoms1,natoms2,-torv,icp.coords);
    }
  }

  double* vx = new double[3];
  vx[1] = vx[2] = 0.;
  vx[0] = 1.0;
  vdw_vector_opt(natoms1,natoms2,vx,icp);
  delete [] vx;

  for (int i=0;i<3*(natoms1+natoms2);i++)
    xyza[i] = icp.coords[i];

  for (int i=0;i<3*natoms1;i++)
    xyz1[i] = xyza[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2[i] = xyza[3*natoms1+i];

  //print_xyz_gen(natoms1+natoms2,icp.anames,xyza);

  int badgeom = 0;
  for (int i=0;i<3*natoms1;i++)
  if (xyz1[i]!=xyz1[i])
  {
    badgeom = 1;
    break;
  }
  for (int i=0;i<3*natoms2;i++)
  if (xyz2[i]!=xyz2[i])
  {
    badgeom = 1;
    break;
  }

  if (badgeom)
  {
    printf(" ERROR: nan detected in add_align, exiting! \n");
    print_xyz_gen(natoms1,anames1,xyz1);
    print_xyz_gen(natoms2,anames2,xyz2);
    exit(1);
  }

  delete [] anames1a;
  delete [] anames2a;
  delete [] xyz1a;
  delete [] xyz2a;

  delete [] v2;
  delete [] v1a;
  delete [] v1b;
  delete [] bonded1;
  delete [] bonded2;
  delete [] bonded3;

  ic1.freemem();
  ic2.freemem();

  delete [] anamesp;
  delete [] anumbersp;
  icp.freemem();


  return;
}

void Align::rotate_around_x(int natoms1, int natoms2, double torv, double* xyz)
{
  //printf("  rotating fragment around x axis: %3.2f (degrees) \n",torv);

  int natoms = natoms1 + natoms2;

  double* xyzn = new double[3*natoms];
  for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;

  double* angles = new double[3];
  double** rotm = new double*[3];
  rotm[0] = new double[3];
  rotm[1] = new double[3];
  rotm[2] = new double[3];

  angles[0] = torv * 3.1415926/180.;
  angles[1] = angles[2] = 0.;

  get_rotation_matrix(rotm,angles);

  for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
  for (int i=0;i<natoms2;i++)
  for (int j=0;j<3;j++)
  for (int k=0;k<3;k++)
    xyzn[3*(natoms1+i)+j] += rotm[j][k]*xyz[3*(i+natoms1)+k];
  for (int i=0;i<3*natoms2;i++)
    xyz[3*natoms1+i] = xyzn[3*natoms1+i];


  delete [] xyzn;
  delete [] angles;
  delete [] rotm[0];
  delete [] rotm[1];
  delete [] rotm[2];
  delete [] rotm;

  return;
}

void Align::align_to_x(int natoms, int t1, int t2, double* xyz, string* anames, int sign, double offset)
{
  //printf(" aligning structure near to x axis, sign: %i offset: %3.2f (degrees) \n",sign,offset);
#if 0
  if (natoms==4)
  {
    printf(" align_to_x, natoms: %i \n",natoms);
    print_xyz_gen(natoms,anames,xyz);
  }
#endif

  //t1 --> 0,0,0
  //t2 --> a,0,0

  offset = offset * 3.14159 / 180.;

  double* x1 = new double[3];
  x1[0] = xyz[3*t2+0] - xyz[3*t1+0];
  x1[1] = xyz[3*t2+1] - xyz[3*t1+1];
  x1[2] = xyz[3*t2+2] - xyz[3*t1+2];

  double n1 = norm(x1,3);
  for (int i=0;i<3;i++) x1[i] = x1[i]/n1;

  //printf(" vector from t1 to t2: %6.4f %6.4f %6.4f \n",x1[0],x1[1],x1[2]);

  double* xyzn = new double[3*natoms];
  for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
  double* xyz1 = new double[3*natoms];
  for (int i=0;i<natoms;i++)
  {
    xyz1[3*i+0] = xyz[3*i+0] - xyz[3*t1+0];
    xyz1[3*i+1] = xyz[3*i+1] - xyz[3*t1+1];
    xyz1[3*i+2] = xyz[3*i+2] - xyz[3*t1+2];
  }
  //print_xyz_gen(natoms,anames,xyz1);
  //printf("b  %s %4.3f %4.3f %4.3f \n",anames[t1].c_str(),xyz1[3*t1+0],xyz1[3*t1+1],xyz1[3*t1+2]);
  //printf("b  %s %4.3f %4.3f %4.3f \n",anames[t2].c_str(),xyz1[3*t2+0],xyz1[3*t2+1],xyz1[3*t2+2]);

  double* angles = new double[3];
  double** rotm = new double*[3];
  rotm[0] = new double[3];
  rotm[1] = new double[3];
  rotm[2] = new double[3];

  double* u1 = new double[3];

 //align second atom to x axis
  x1[0] = xyz1[3*t2+0] - xyz1[3*t1+0];
  x1[1] = xyz1[3*t2+1] - xyz1[3*t1+1];
  x1[2] = xyz1[3*t2+2] - xyz1[3*t1+2];

//  n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
  n1 = sqrt(x1[0]*x1[0]+x1[2]*x1[2]);
  if (n1>0.)
  for (int i=0;i<3;i++) x1[i] = x1[i]/n1;

  u1[0] = x1[0];
  u1[1] = x1[1];
  u1[2] = x1[2];
  //printf(" u1: %4.3f %4.3f %4.3f \n",u1[0],u1[1],u1[2]);

  angles[0] = angles[1] = angles[2] = 0.;
  angles[1] = acos(u1[0])+offset/3.14159;
  if (u1[2]<0.) angles[1] = angles[1] * -1;
  //angles[2] = acos(u1[1]);
  //printf(" angle of %i %i to x axis: %4.3f \n",t1,t2,angles[1]/3.14159);

  get_rotation_matrix(rotm,angles);

  for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
  for (int i=0;i<natoms;i++)
  for (int j=0;j<3;j++)
  for (int k=0;k<3;k++)
    xyzn[3*i+j] += rotm[j][k]*xyz1[3*i+k];
  for (int i=0;i<3*natoms;i++)
    xyz1[i] = xyzn[i];

  //print_xyz_gen(natoms,anames,xyzn);
  //printf("1  %s %4.3f %4.3f %4.3f \n",anames[t1].c_str(),xyz1[3*t1+0],xyz1[3*t1+1],xyz1[3*t1+2]);
  //printf("1  %s %4.3f %4.3f %4.3f \n",anames[t2].c_str(),xyz1[3*t2+0],xyz1[3*t2+1],xyz1[3*t2+2]);


 //start second rotation
  x1[0] = xyz1[3*t2+0] - xyz1[3*t1+0];
  x1[1] = xyz1[3*t2+1] - xyz1[3*t1+1];
  x1[2] = xyz1[3*t2+2] - xyz1[3*t1+2];

//  n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
  n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]);
  if (n1>0.0)
  for (int i=0;i<3;i++) x1[i] = x1[i]/n1;

  u1[0] = x1[0];
  u1[1] = x1[1];
  u1[2] = x1[2];
  //printf(" u1: %4.3f %4.3f %4.3f \n",u1[0],u1[1],u1[2]);

  angles[0] = angles[1] = angles[2] = 0.;
  angles[2] = acos(u1[0])+offset;
  if (u1[1]>0.) angles[2] = angles[2] * -1;
  //printf(" angle of %i %i to x axis: %4.3f \n",t1,t2,angles[2]/3.14159);

  get_rotation_matrix(rotm,angles);

  for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
  for (int i=0;i<natoms;i++)
  for (int j=0;j<3;j++)
  for (int k=0;k<3;k++)
    xyzn[3*i+j] += rotm[j][k]*xyz1[3*i+k];
  for (int i=0;i<3*natoms;i++)
    xyz1[i] = xyzn[i];

//  if (natoms==4)
//    print_xyz_gen(natoms,anames,xyzn);
  //printf("2  %s %4.3f %4.3f %4.3f \n",anames[t1].c_str(),xyz1[3*t1+0],xyz1[3*t1+1],xyz1[3*t1+2]);
  //printf("2  %s %4.3f %4.3f %4.3f \n",anames[t2].c_str(),xyz1[3*t2+0],xyz1[3*t2+1],xyz1[3*t2+2]);

  if (sign==-1)
  {
    angles[0] = angles[1] = angles[2] = 0.;
    angles[2] = 3.14159265;
    //if (u1[1]>0.) angles[2] = angles[2] * -1;
    //printf(" angle of %i %i to x axis: %4.3f \n",t1,t2,angles[2]/3.14159);

    get_rotation_matrix(rotm,angles);

    for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<3;j++)
    for (int k=0;k<3;k++)
      xyzn[3*i+j] += rotm[j][k]*xyz1[3*i+k];
    for (int i=0;i<3*natoms;i++)
      xyz1[i] = xyzn[i];
  }
  

#if 0
  //checking final geom
  double THRESH = 0.1;
  if (abs(xyzn[3*t1+0]) > THRESH 
   || abs(xyzn[3*t1+1]) > THRESH
   || abs(xyzn[3*t1+2]) > THRESH
   || abs(xyzn[3*t2+1]) > THRESH
   || abs(xyzn[3*t2+2]) > THRESH)
  {
    printf("\n WARNING: align did not set geom to 0,0,0/a,0,0 \n");
    //align_to_x(natoms,t1,t2,xyzn,anames);
  }
#endif

  for (int i=0;i<3*natoms;i++)
    xyz[i] = xyzn[i];

  delete [] angles;
  for (int i=0;i<3;i++)
    delete [] rotm[i];
  delete [] rotm;

  delete [] xyzn;
  delete [] xyz1;

  return;
}


void Align::align_v1(int nvf, double* v1)
{
  double dot;
  for (int i=1;i<nvf;i++)
  {
    dot = 0.;
    for (int j=0;j<3;j++)
      dot += v1[j] * v1[3*i+j];
    //printf("  dot (%i to 0): %4.3f \n",i,dot);
    if (dot<0.0)
    for (int j=0;j<3;j++)
      v1[3*i+j] = - v1[3*i+j];
  }

  return;
}

void Align::linear_right(double* v1, int a1, int* bonded, double* xyz)
{
  double* x1 = new double[3];
  for (int i=0;i<3;i++)
    x1[i] = xyz[3*a1+i] - xyz[3*bonded[0]+i];

  double* x2 = new double[3];
  x2[0] = x2[1] = 0.;
  x2[2] = 1.;
 
  cross(v1,x1,x2);

  double n1 = norm(v1,3);
  if (n1<0.001)
  {
    x2[0] = x2[2] = 0.;
    x2[1] = 1.0;
    cross(v1,x1,x2);
  }

  delete [] x1;
  delete [] x2;

  return;
}

void Align::planar_cross(double* v1, int a1, int* bonded, double* xyz)
{
  double* x1 = new double[3];
  double* x2 = new double[3];

  int b2 = bonded[0];
  int b3 = bonded[1];
  for (int i=0;i<3;i++)
    x1[i] = xyz[3*b2+i] - xyz[3*a1+i];
  for (int i=0;i<3;i++)
    x2[i] = xyz[3*b3+i] - xyz[3*a1+i];

  double n1 = norm(x1,3);
  for (int i=0;i<3;i++)
    x1[i] = x1[i] / n1;
  n1 = norm(x2,3);
  for (int i=0;i<3;i++)
    x2[i] = x2[i] / n1;

  cross(v1,x1,x2);

  delete [] x1;
  delete [] x2;


  return;
}


int Align::get_bonds(int a1, ICoord ic1, int* bonded)
{
  int nfound = 0;

  for (int i=0;i<ic1.natoms;i++)
  if (i!=a1 && ic1.bond_exists(i,a1))
    bonded[nfound++] = i;

#if 0
  for (int i=0;i<nfound;i++)
    printf(" %2i bonded to %2i \n",a1+1,bonded[i]+1);
#endif

  return nfound;
}

int Align::check_frag(int a1, int a2)
{
  int f1,f2;
  if (a1<natoms1) f1 = 1;
  else f1 = 2;
  if (a2<natoms1) f2 = 1;
  else f2 = 2;

  int same_frag = 1;
  if (f1!=f2) same_frag = 0;

  return same_frag;
}

int Align::check_frag_3(int a1, int a2)
{
  int f1,f2;
  if (a1<natoms1+natoms2) f1 = 1;
  else f1 = 2;
  if (a2<natoms1+natoms2) f2 = 1;
  else f2 = 2;

  int same_frag = 1;
  if (f1!=f2) same_frag = 0;

  return same_frag;
}


void Align::init(int natoms1i, string* anames1i, int* anumbers1i, double* xyz1i, int natoms2i, string* anames2i, int* anumbers2i, double* xyz2i)
{
  //printf(" Align init: %i %i \n",natoms1i, natoms2i);
  inited = 1;

  natoms1 = natoms1i;
  natoms2 = natoms2i;
  anames1 = new string[natoms1];
  anames2 = new string[natoms2];
  anumbers1 = new int[natoms1];
  anumbers2 = new int[natoms2];
  xyz1 = new double[3*natoms1];
  xyz2 = new double[3*natoms2];
  xyz1s = new double[3*natoms1];
  xyz2s = new double[3*natoms2];
  xyza = new double[3*(natoms1+natoms2)];

  for (int i=0;i<natoms1;i++)
  {
    anames1[i] = anames1i[i];
    anumbers1[i] = anumbers1i[i];
  }
  for (int i=0;i<natoms2;i++)
  {
    anames2[i] = anames2i[i];
    anumbers2[i] = anumbers2i[i];
  }
  for (int i=0;i<3*natoms1;i++)
    xyz1s[i] = xyz1[i] = xyz1i[i];
  for (int i=0;i<3*natoms2;i++)
    xyz2s[i] = xyz2[i] = xyz2i[i];
  for (int i=0;i<3*natoms1;i++)
    xyza[i] = xyz1i[i];
  for (int i=0;i<3*natoms2;i++)
    xyza[3*natoms1+i] = xyz2i[i];

  anames3 = NULL;
  anumbers3 = NULL;
  xyz3 = NULL;
  xyza3 = NULL;

  return;
}

void Align::freemem()
{
  //printf(" Align free: %i %i \n",natoms1,natoms2);

  if (inited==1 || inited==3)
  {
    if (anames1!=NULL) delete [] anames1;
    if (anames2!=NULL) delete [] anames2;
    if (anumbers1!=NULL) delete [] anumbers1;
    if (anumbers2!=NULL) delete [] anumbers2;
    if (xyz1!=NULL) delete [] xyz1;
    if (xyz2!=NULL) delete [] xyz2;
    if (xyza!=NULL) delete [] xyza;
  }
  if (inited==3)
  {
    if (anames3!=NULL) delete [] anames3;
    if (anumbers3!=NULL) delete [] anumbers3;
    if (xyz3!=NULL) delete [] xyz3;
    if (xyza3!=NULL) delete [] xyza3;
  }

  inited = 0;

  return;
}

void Align::print_xyz()
{
  printf(" %2i \n\n",natoms1+natoms2);
  for (int i=0;i<natoms1;i++)
    printf(" %s %9.6f %9.6f %9.6f \n",anames1[i].c_str(),xyz1[3*i+0],xyz1[3*i+1],xyz1[3*i+2]);
  for (int i=0;i<natoms2;i++)
    printf(" %s %9.6f %9.6f %9.6f \n",anames2[i].c_str(),xyz2[3*i+0],xyz2[3*i+1],xyz2[3*i+2]);
  printf("\n");

  return;
}

void Align::print_xyz_3()
{
  if (inited!=3)
  {
    printf(" WARNING: inited is not 3 in print_xyz_3 \n");
    return;
  }

  printf(" %2i \n\n",natoms1+natoms2+natoms3);
  for (int i=0;i<natoms1;i++)
    printf(" %s %9.6f %9.6f %9.6f \n",anames1[i].c_str(),xyz1[3*i+0],xyz1[3*i+1],xyz1[3*i+2]);
  for (int i=0;i<natoms2;i++)
    printf(" %s %9.6f %9.6f %9.6f \n",anames2[i].c_str(),xyz2[3*i+0],xyz2[3*i+1],xyz2[3*i+2]);
  for (int i=0;i<natoms3;i++)
    printf(" %s %9.6f %9.6f %9.6f \n",anames3[i].c_str(),xyz3[3*i+0],xyz3[3*i+1],xyz3[3*i+2]);
  printf("\n");

  return;
}
