// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "zstruct.h"
#include "utils.h"
#include "constants.h"
#include "print.h"



int ZStruct::diff_structureigq(int natoms1, int natoms2, string* anames1, string* anames2, int* anumbers1, int* anumbers2, double* xyz1, double* xyz2)
{
  if (natoms1!=natoms2) return 1;

  int diff = 0;
  int cont = 1;
  int natoms0 = natoms1;

#if 0
  printf(" in diff_structureigq \n");
  print_xyz_gen(natoms0,anames1,xyz1);
  print_xyz_gen(natoms0,anames2,xyz2);
#endif

 //make sure both structures have same # of each element
  int* anumc1 = new int[nelem];
  int* anumc2 = new int[nelem];
  for (int i=0;i<nelem;i++)
    anumc1[i] = anumc2[i] = 0;
  for (int i=0;i<natoms0;i++)
    anumc1[anumbers1[i]]++;
  for (int i=0;i<natoms0;i++)
    anumc2[anumbers2[i]]++;
  for (int i=1;i<nelem;i++)
  if (anumc1[i]!=anumc2[i])
    diff++;
  if (diff>0) cont = 0;

  if (!cont)
  {
    delete [] anumc1;
    delete [] anumc2;
    return diff;
  }

#if 0
  for (int i=0;i<nelem;i++)
  if (anumc1[i])
    printf(" element: %i \n",i);
#endif
 

 //arrange elements/xyz's into standard order
  double* xyz1a = new double[3*natoms0];
  double* xyz2a = new double[3*natoms0];
  string* anames0 = new string[natoms0];
  int* anumbers0 = new int[natoms0];

  int naf = 0;
  for (int i=1;i<nelem;i++)
  if (anumc1[i])
  {
    for (int j=0;j<natoms0;j++)
    if (anumbers1[j]==i)
    {
      xyz1a[3*naf+0] = xyz1[3*j+0];
      xyz1a[3*naf+1] = xyz1[3*j+1];
      xyz1a[3*naf+2] = xyz1[3*j+2];
      anames0[naf]    = anames1[j];
      anumbers0[naf]  = anumbers1[j];
      naf++;
    }
  }
  naf = 0;
  for (int i=1;i<nelem;i++)
  if (anumc1[i])
  {
    for (int j=0;j<natoms0;j++)
    if (anumbers2[j]==i)
    {
      xyz2a[3*naf+0] = xyz2[3*j+0];
      xyz2a[3*naf+1] = xyz2[3*j+1];
      xyz2a[3*naf+2] = xyz2[3*j+2];
      naf++;
    }
  }

#if 0
  printf(" printing xyz/a after arrange \n");
  printf(" %i \n\n",natoms0);
  for (int i=0;i<natoms0;i++)
    printf(" %s %6.5f %6.5f %6.5f \n",anames0[i].c_str(),xyz1a[3*i+0],xyz1a[3*i+1],xyz1a[3*i+2]);
  printf(" %i \n\n",natoms0);
  for (int i=0;i<natoms0;i++)
    printf(" %s %6.5f %6.5f %6.5f \n",anames0[i].c_str(),xyz2a[3*i+0],xyz2a[3*i+1],xyz2a[3*i+2]);
#endif

  diff = diff_structureiq(natoms0,anames0,anumbers0,xyz1a,xyz2a);

  delete [] anumc1;
  delete [] anumc2;

  delete [] xyz1a;
  delete [] xyz2a;
  delete [] anames0;
  delete [] anumbers0;

  return diff;
}

int ZStruct::diff_structureiq(int natoms1, string* anames1, int* anumbers1, double* xyz1, double* xyz2)
{
  //assumes both xyz's have same # of atoms
  //assumes same element ordering

#if 0
  //debug
  if (xyz1==NULL || xyz2==NULL)
    printf("  xyz problem in diff_structureiq \n");
  if (anames1==NULL || anumbers1==NULL)
    printf("  anames or anumbers problem in diff_structureiq \n");
#endif

  //printf("\n\n diff_structureiq test \n");

  int diff = 0;
  double* xyz3 = new double[3*natoms1];
  for (int i=0;i<3*natoms1;i++)
    xyz3[i] = xyz2[i];

  ICoord test1;
  ICoord test2;
  test1.alloc(natoms1);
  test2.alloc(natoms1);
  test1.farBond = 1.1;
  test2.farBond = 1.1;
  test1.reset(natoms1,anames1,anumbers1,xyz1);
  test2.reset(natoms1,anames1,anumbers1,xyz2);
  test1.ic_create();
  test2.ic_create();

  int cont = 1;

// diff struct without rotation
  if (diff_structureq(natoms1,anames1,anumbers1,xyz1,xyz2))
  {
    //printf(" first test no info \n");
  }
  else
  {
    //printf(" same found \n");
    diff = 0;
    cont = 0;
  }

// make element list
  int* atype = new int[natoms1];
  for (int i=0;i<natoms1;i++) 
    atype[i] = -1;
  atype[0] = anumbers1[0];
  int ntype = 1;
  for (int i=0;i<natoms1;i++) 
  {
    int found = 0;
    for (int j=0;j<ntype;j++) 
    if (anumbers1[i]==atype[j])
      found = 1;
    if (!found)
      atype[ntype++] = anumbers1[i];
  }  
#if 0
  printf(" printing element types:");
  for (int i=0;i<ntype;i++)
    printf(" %i",atype[i]);
  printf("\n");
#endif

// count number of each coordn for each element
  int* coordn1 = new int[12];
  int* coordn2 = new int[12];
  for (int i=0;i<ntype;i++)
  {
    for (int j=0;j<12;j++) coordn1[j] = 0;
    for (int j=0;j<12;j++) coordn2[j] = 0;
    for (int j=0;j<natoms1;j++)
    if (anumbers1[j] == atype[i])
    {
      coordn1[test1.coordn[j]]++;
      coordn2[test2.coordn[j]]++;
    }
//    for (int j=0;j<12;j++)
//      printf(" AtNum: %i coord#: %i coordn1: %i coordn2: %i \n",atype[i],j,coordn2[j],coordn2[j]);
    for (int j=0;j<12;j++)
    if (coordn1[j]!=coordn2[j])
      diff++;
  }
  //printf(" after checking elements for same coordn, diff: %i \n",diff);
  if (diff>0) cont = 0;

//  if (cont)
//  printf(" diff_structureiq element swapping \n");
//  printf("h1"); fflush(stdout);

// rotate heavy atoms by coordn
  if (cont)
  for (int i=0;i<ntype;i++)
  if (atype[i]!=1)
  {
    for (int j=0;j<natoms1;j++)
    if (anumbers1[j]==atype[i])
    if (test1.coordn[j]!=test2.coordn[j])
    {
      //printf(" need a swap partner for %i \n",j);
      int s1 = -1;
      int cn = test1.coordn[j];
      for (int k=0;k<natoms1;k++)
      if (j!=k && cn==test2.coordn[k] && anumbers1[k]==atype[i])
      {
        s1 = k;
        //printf(" found: %i \n",s1);
        break;
      }
      swap_atoms(test2.coords,j,s1);
      test2.ic_create();
    }
  }

  int* clist1 = new int[6*natoms1];
  int* clist2 = new int[6*natoms1];
  for (int i=0;i<6*natoms1;i++) clist1[i] = -1;
  for (int i=0;i<6*natoms1;i++) clist2[i] = -1;
  if (cont)
  for (int i=0;i<natoms1;i++)
  {
    int nconn1 = 0;
    int nconn2 = 0;
    for (int j=0;j<test1.nbonds;j++)
    {
      if (test1.bonds[j][0] == i)
        clist1[6*i+nconn1++] = anumbers1[test1.bonds[j][1]];
      else if (test1.bonds[j][1] == i)
        clist1[6*i+nconn1++] = anumbers1[test1.bonds[j][0]];
      if (test2.bonds[j][0] == i)
        clist2[6*i+nconn2++] = anumbers1[test2.bonds[j][1]];
      else if (test2.bonds[j][1] == i)
        clist2[6*i+nconn2++] = anumbers1[test2.bonds[j][0]];
    }
  }
  if (cont)
  for (int i=0;i<natoms1;i++)
  {
    for (int j=0;j<5;j++)
    if (clist1[6*i+j]<clist1[6*i+j+1])
    {
      //printf("  if %i %i %i %i \n",i,j,clist1[6*i+j],clist1[6*i+j+1]);
      int tmp = clist1[6*i+j];
      clist1[6*i+j] = clist1[6*i+j+1];
      clist1[6*i+j+1] = tmp;
      j = -1;
    }
    for (int j=0;j<5;j++)
    if (clist2[6*i+j]<clist2[6*i+j+1])
    {
      //printf("  if %i %i %i %i \n",i,j,clist2[6*i+j],clist2[6*i+j+1]);
      int tmp = clist2[6*i+j];
      clist2[6*i+j] = clist2[6*i+j+1];
      clist2[6*i+j+1] = tmp;
      j = -1;
    }
  }
#if 0
  for (int i=0;i<natoms1;i++)
  {
    printf(" atom %i has element connections: %i %i %i %i %i %i \n",i,clist1[6*i+0],clist1[6*i+1],clist1[6*i+2],clist1[6*i+3],clist1[6*i+4],clist1[6*i+5]);
    printf(" atom %i has element connections: %i %i %i %i %i %i \n",i,clist2[6*i+0],clist2[6*i+1],clist2[6*i+2],clist2[6*i+3],clist2[6*i+4],clist2[6*i+5]);
  }
#endif

  //printf("h2"); fflush(stdout);

  if (diff>0) cont = 0;
  int* fixed = new int[natoms1];
  for (int i=0;i<natoms1;i++)
    fixed[i] = test1.coordn[i];

  int* bondok = new int[test1.nbonds];
  for (int i=0;i<test1.nbonds;i++) bondok[i] = 0;

// look for same bonds
  int nsamecoordn = 0;
  int found = 0;
  if (cont)
  for (int i=0;i<test1.nbonds;i++)
  {
    if (test2.bond_exists(test1.bonds[i][0],test1.bonds[i][1]))
    {
//      printf(" found bond: %i %i \n",test1.bonds[i][0],test1.bonds[i][1]);
      fixed[test1.bonds[i][0]]--;
      fixed[test1.bonds[i][1]]--;
      bondok[i] = 1;
    }
  }

  //printf("h3"); fflush(stdout);

// heavy atoms not fully matched to bonds, check local attachment
  int nsame = 0;
  if (cont)
  for (int i=0;i<ntype;i++)
  if (atype[i]!=1)
  {
    for (int j=0;j<natoms1;j++)
    if (anumbers1[j]==atype[i] && fixed[j]>0)
    {
      int same1 = 1; 
      for (int l=0;l<6;l++)
      if (clist1[6*j+l]!=clist2[6*j+l])
        same1 = 0;
      if (same1)
      {
        //printf(" atom %i same elemental coord sphere \n",j);
        fixed[j] = 0;
      }
      if (!same1)
      {
        int found = 0;
        for (int k=0;k<natoms1;k++)
        if (j!=k && anumbers1[k]==atype[i] && fixed[k]>0)
        {
          int same2 = 1;
          for (int l=0;l<6;l++)
          if (clist1[6*j+l]!=clist2[6*k+l])
            same2 = 0;
          if (same2)
          {
            found = 1;
            //printf(" atoms %i %i could be swapped \n",j,k);
            fixed[k] = 0;
            break;
          }
        }
        if (!found) //CPMZ heavy atom elemental coord sphere criterion 
          diff++;
      } //if !same1
    } //loop over heavy atom j
  } //loop over elements i

  if (diff>0) cont = 0;
#if 0
  if (cont)
  for (int i=0;i<natoms1;i++)
    printf(" test1.coordn[%i]: %i fixed[%i]: %i \n",i,test1.coordn[i],i,fixed[i]);
#endif

  if (cont)
  {
   // test1.print_bonds();
   // test2.print_bonds();
   // test1.print_xyz();
   // test2.print_xyz();
  }

#if 0
//not using?
  int* ne = new int[ntype];
  for (int i=0;i<ntype;i++)
    ne[i] = 0;
  int** elemlist = new int*[ntype];
  for (int i=0;i<ntype;i++)
    elemlist[i] = new int[natoms1];
  for (int i=0;i<ntype;i++)
  for (int j=0;j<natoms1;j++)
  if (anumbers1[j] == atype[i])
    elemlist[i][ne[i]++] = j;
#endif


//rotate within same element then compare...
#if 0
//not fully implemented
  if (cont)
  for (int b=0;b<test1.nbonds;b++)
  if (!bondok[b])
  {
    int a1 = test1.bonds[b][0];
    int a2 = test1.bonds[b][1];
   // printf(" looking for bond %i(%s) %i(%s) \n",a1,anames[a1].c_str(),a2,anames[a2].c_str());
   // printf(" fixed[a1]: %i fixed[a2]: %i \n",fixed[a1],fixed[a2]);
    int a12 = -1;
    int a21 = -1;
    if (fixed[a1]>0 && anumbers1[a1]==1)
    {
      a12 = a1;
      a21 = a2;
    }
    else if (fixed[a2]>0 && anumbers1[a2]==1)
    {
      a12 = a2;
      a21 = a1;
    }
    if (a12!=-1)
    for (int i=0;i<natoms1;i++)
    if (i!=a12 && fixed[i]>0 && anumbers1[i]==anumbers1[a12])
    {
      //printf(" swap candidate: %i %i \n",a12,i);
      if (test2.bond_exists(i,a21))
      {
        printf(" swap okay: %i %i becomes %i %i \n",a12,a21,i,a21);
//        fixed[i]--;
//        fixed[a21]--;
        bondok[b] = 1;
        break;
      }
    } //if candidate i for a12
  } //loop over not okay bonds

  if (cont)
  for (int i=0;i<test1.nbonds;i++)
  if (!bondok[i])
  {
    //printf(" bond %i not found \n",i);
    diff++;
  }
  printf(" printing which atoms may be rotated \n");
  for (int i=0;i<natoms1;i++)
    printf(" test1.coordn[%i]: %i fixed[%i]: %i \n",i,test1.coordn[i],i,fixed[i]);
#endif

  //printf("hd"); fflush(stdout);

  delete [] atype;
  delete [] clist1;
  delete [] clist2;
  delete [] fixed;
  delete [] bondok;
  delete [] coordn1;
  delete [] coordn2;
  delete [] xyz3;
  test1.freemem();
  test2.freemem();
  
  return diff;
}

void ZStruct::swap_atoms(double* xyz, int a1, int a2) {

  double* xyz3 = new double[3];

  xyz3[0] = xyz[3*a1+0];
  xyz3[1] = xyz[3*a1+1];
  xyz3[2] = xyz[3*a1+2];
  xyz[3*a1+0] = xyz[3*a2+0];
  xyz[3*a1+1] = xyz[3*a2+1];
  xyz[3*a1+2] = xyz[3*a2+2];
  xyz[3*a2+0] = xyz3[0];
  xyz[3*a2+1] = xyz3[1];
  xyz[3*a2+2] = xyz3[2];

  delete [] xyz3;
  return;
}

//uses coordination # only
int ZStruct::diff_structurec(int natoms1, string* anames1, int* anumbers1, double* xyz1, double* xyz2)
{
  int diff = 0;

  ICoord test1;
  ICoord test2;
  test1.alloc(natoms1);
  test2.alloc(natoms1);
  test1.farBond = 1.1;
  test2.farBond = 1.1;
  test1.reset(natoms1,anames1,anumbers1,xyz1);
  test2.reset(natoms1,anames1,anumbers1,xyz2);
  test1.ic_create();
  test2.ic_create();

  //printf(" testing diff_structurecq \n");

  int nsamecoordn = 0;
  for (int k=0;k<natoms1;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms1) 
  {
    printf(" same coordination #s \n");
    diff = 0;
  }
  else
  {
    diff = natoms1-nsamecoordn;
    printf(" different coordination# (by %i) \n",diff);
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}


//uses coordination # only
int ZStruct::diff_structurecq(int natoms1, string* anames1, int* anumbers1, double* xyz1, double* xyz2)
{
  int diff = 0;
  ICoord test1;
  ICoord test2;
  test1.alloc(natoms1);
  test2.alloc(natoms1);
  test1.farBond = 1.1;
  test2.farBond = 1.1;
  test1.reset(natoms1,anames1,anumbers1,xyz1);
  test2.reset(natoms1,anames1,anumbers1,xyz2);
  test1.ic_create();
  test2.ic_create();

  //printf(" testing diff_structurecq \n");

  int nsamecoordn = 0;
  for (int k=0;k<natoms1;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms1) 
  {
    //printf(" same coordination #s \n");
    diff = 0;
  }
  else
  {
    diff = natoms1-nsamecoordn;
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}



int ZStruct::diff_structure(int natoms1, string* anames1, int* anumbers1, double* xyz1, double* xyz2)
{
  //printf(" in diff_structure for natoms1: %i \n",natoms1); fflush(stdout);

  int diff = 0;

  ICoord test1;
  ICoord test2;
  test1.alloc(natoms1);
  test2.alloc(natoms1);
  test1.farBond = 1.1;
  test2.farBond = 1.1;
  test1.reset(natoms1,anames1,anumbers1,xyz1);
  test2.reset(natoms1,anames1,anumbers1,xyz2);
  test1.ic_create();
  test2.ic_create();

  int nsamecoordn = 0;
  for (int k=0;k<natoms1;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms1) 
  {
    //printf(" same coordination #s, checking connectivity \n");
    int found = 0;
    for (int i=0;i<test1.nbonds;i++)
    {
      found = 0;
      if (test2.bond_exists(test1.bonds[i][0],test1.bonds[i][1]))
        found = 1;
      else
        break;
    }
    if (found) diff = 0;
    else diff = 2;
  }
  else
  {
    diff = natoms1-nsamecoordn;
//    printf(" different coordination# (by %i) \n",diff);
//    nmopacdiff++;
//    test1.print_xyz();
//    test2.print_xyz();
//    test1.print_bonds();
//    test2.print_bonds();
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}


int ZStruct::diff_structureq(int natoms1, string* anames1, int* anumbers1, double* xyz1, double* xyz2)
{ 
  int diff = 0;
  ICoord test1;
  ICoord test2;
  test1.alloc(natoms1);
  test2.alloc(natoms1);
  test1.farBond = 1.1;
  test2.farBond = 1.1;
  test1.reset(natoms1,anames1,anumbers1,xyz1);
  test2.reset(natoms1,anames1,anumbers1,xyz2);
  test1.ic_create();
  test2.ic_create();

  //printf(" testing diff_structure \n");

  int nsamecoordn = 0;
  for (int k=0;k<natoms1;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms1) 
  {
    //printf(" same coordination #s, checking connectivity \n");
    int found = 0;
    for (int i=0;i<test1.nbonds;i++)
    {
      found = 0;
      if (test2.bond_exists(test1.bonds[i][0],test1.bonds[i][1]))
        found = 1;
      else
        break;
    }
    if (found) diff = 0;
    else diff = 2;
  }
  else
  {
    diff = natoms1-nsamecoordn;
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}

