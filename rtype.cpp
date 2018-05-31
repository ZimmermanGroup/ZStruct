// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "rtype.h"


void RTYPE::print()
{
  printf("   natoms/nadd/nbrk: %i %i %i   Pr: %4.3f id: %4i \n",natoms,nadd,nbrks,Pr,id);
  printf("    adds:");
  for (int i=0;i<nadd;i++)
  if (coordn[add[2*i+0]]>-1 || coordn[add[2*i+1]]>-1)
    printf("  %s%i - %s%i",anames[add[2*i+0]].c_str(),coordn[add[2*i+0]],anames[add[2*i+1]].c_str(),coordn[add[2*i+1]]);
  else
    printf("  %s - %s",anames[add[2*i+0]].c_str(),anames[add[2*i+1]].c_str());
  printf("\n");
  printf("    brks:");
  for (int i=0;i<nbrks;i++)
  if (coordn[brks[2*i+0]]>-1 || coordn[brks[2*i+1]]>-1)
    printf("  %s%i - %s%i",anames[brks[2*i+0]].c_str(),coordn[brks[2*i+0]],anames[brks[2*i+1]].c_str(),coordn[brks[2*i+1]]);
  else
    printf("  %s - %s",anames[brks[2*i+0]].c_str(),anames[brks[2*i+1]].c_str());
  printf("\n");


  return;
}

int RTYPE::match(int natoms1, int* anumbers1, int* coordn1, int nadd1, int* add1, int nbrks1, int* brks1)
{
  if (natoms1!=natoms)
    return 0;

  if (nadd1!=nadd || nbrks1!=nbrks)
    return 0;


  sort_rxn1(nadd1,add1,nbrks1,brks1,anumbers1,coordn1);

  //printf("   in match: nadd1/nbrks1: %i %i  (id: %4i) \n",nadd1,nbrks1,id);

  int* ab_used = new int[nadd+nbrks];
  int nfound = 0;

  for (int i=0;i<nadd+nbrks;i++) ab_used[i] = 0;
  for (int i=0;i<nadd1;i++)
  if (nfound>-1)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int e1 = anumbers1[a1];
    int e2 = anumbers1[a2];
    int c1 = coordn1[a1];
    int c2 = coordn1[a2];

    //printf("    add comparing: e1/e2: %i %i \n",e1,e2);
    int found = 0;
    for (int j=0;j<nadd;j++)
    if (ab_used[j]==0)
//    if ((e1==anumbers[add[2*j+0]] || anumbers[add[2*j+0]]==-1) && (e2==anumbers[add[2*j+1]] || anumbers[add[2*j+1]]==-1))
//    if ((c1==coordn[add[2*j+0]] || coordn[add[2*j+0]]==-1) && (c2==coordn[add[2*j+1]] || coordn[add[2*j+1]]==-1))
    if ((e1==anumbers[add[2*j+0]] || e1==-1 || anumbers[add[2*j+0]]==-1) && (e2==anumbers[add[2*j+1]] || e2==-1 || anumbers[add[2*j+1]]==-1))
    if ((c1==coordn[add[2*j+0]] || c1==-1 || coordn[add[2*j+0]]==-1) && (c2==coordn[add[2*j+1]] || c2==-1 || coordn[add[2*j+1]]==-1))
    {
      //printf("    add matched: e1/e2: %i %i \n",e1,e2);
      ab_used[j] = 1;
      found = 1;
      nfound++;
      break;
    }
    if (!found)
    {
      //printf("    not found: e1/e2: %i %i \n",e1,e2);
      nfound = -1;
      break;
    }
  } //loop i over nadd1

  for (int i=0;i<nadd+nbrks;i++) ab_used[i] = 0;
  for (int i=0;i<nbrks1;i++)
  if (nfound>-1)
  {
    int b1 = brks1[2*i+0];
    int b2 = brks1[2*i+1];
    int e1 = anumbers1[b1];
    int e2 = anumbers1[b2];
    int c1 = coordn1[b1];
    int c2 = coordn1[b2];

    //printf("    brk comparing: e1/e2: %i %i \n",e1,e2);
    int found = 0;
    for (int j=0;j<nbrks;j++)
    if (ab_used[j]==0)
//    if ((e1==anumbers[brks[2*j+0]] || anumbers[brks[2*j+0]]==-1) && (e2==anumbers[brks[2*j+1]] || anumbers[brks[2*j+1]]==-1))
//    if ((c1==coordn[brks[2*j+0]] || coordn[brks[2*j+0]]==-1) && (c2==coordn[brks[2*j+1]] || coordn[brks[2*j+1]]==-1))
    if ((e1==anumbers[brks[2*j+0]] || e1==-1 || anumbers[brks[2*j+0]]==-1) && (e2==anumbers[brks[2*j+1]] || e2==-1 || anumbers[brks[2*j+1]]==-1))
    if ((c1==coordn[brks[2*j+0]] || c1==-1 || coordn[brks[2*j+0]]==-1) && (c2==coordn[brks[2*j+1]] || c2==-1 || coordn[brks[2*j+1]]==-1))
    {
      //printf("    brk matched: e1/e2: %i %i \n",e1,e2);
      ab_used[j] = 1;
      found = 1;
      nfound++;
      break;
    }
    if (!found)
    {
      //printf("    not found: e1/e2: %i %i \n",e1,e2);
      nfound = -1;
      break;
    }
  } //loop i over nadd1


  delete [] ab_used;


  int success = 1;
  if (nfound<0) success = 0;

  return success;
}



void RTYPE::sort_rxn1(int nadd1, int* add1, int nbrks1, int* brks1, int* anumbers1, int* coordn1)
{
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int e1 = anumbers1[a1];
    int e2 = anumbers1[a2];
    int c1 = coordn1[a1];
    int c2 = coordn1[a2];
    if (e2>e1)
    {
      add1[2*i+0] = a2;
      add1[2*i+1] = a1;
    }
    else if (e1==e2 && c2>c1)
    {
      add1[2*i+0] = a2;
      add1[2*i+1] = a1;      
    }
  }
  for (int i=0;i<nbrks1;i++)
  {
    int b1 = brks1[2*i+0];
    int b2 = brks1[2*i+1];
    int e1 = anumbers1[b1];
    int e2 = anumbers1[b2];
    int c1 = coordn1[b1];
    int c2 = coordn1[b2];

    if (e2>e1)
    {
      brks1[2*i+0] = b2;
      brks1[2*i+1] = b1;
    }
    else if (e1==e2 && c2>c1)
    {
      brks1[2*i+0] = b2;
      brks1[2*i+1] = b1;
    }
  }
 
  return;
}

void RTYPE::sort_rxn()
{
  for (int i=0;i<nadd;i++)
  {
    int a1 = add[2*i+0];
    int a2 = add[2*i+1];
    int e1 = anumbers[a1];
    int e2 = anumbers[a2];
    int c1 = coordn[a1];
    int c2 = coordn[a2];

    if (e2>e1)
    {
      add[2*i+0] = a2;
      add[2*i+1] = a1;
    }
    else if (e1==e2 && c2>c1)
    {
      add[2*i+0] = a2;
      add[2*i+1] = a1;
    }
  }
  for (int i=0;i<nbrks;i++)
  {
    int b1 = brks[2*i+0];
    int b2 = brks[2*i+1];
    int e1 = anumbers[b1];
    int e2 = anumbers[b2];
    int c1 = coordn[b1];
    int c2 = coordn[b2];

    if (e2>e1)
    {
      brks[2*i+0] = b2;
      brks[2*i+1] = b1;
    }
    else if (e1==e2 && c2>c1)
    {
      brks[2*i+0] = b2;
      brks[2*i+1] = b1;
    }
  }

  return;
}

void RTYPE::set_rxn(int natoms1, string* anames1, int* anumbers1, int* coordn1, int nadd1, int* add1, int nbrks1, int* brks1)
{
  if (alloced==1) freemem();

  nadd = nadd1;
  nbrks = nbrks1;
  add = new int[2*nadd];
  brks = new int[2*nbrks];
  for (int i=0;i<2*nadd;i++)
    add[i] = add1[i];
  for (int i=0;i<2*nbrks;i++)
    brks[i] = brks1[i];

  natoms = natoms1;
  anames = new string[natoms];
  anumbers = new int[natoms];
  coordn = new int[natoms];
  for (int i=0;i<natoms;i++)
    anames[i] = anames1[i];
  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers1[i];
  for (int i=0;i<natoms;i++)
    coordn[i] = coordn1[i];

  sort_rxn();

  return;
}


double RTYPE::get_pr()
{
  return Pr;
}

void RTYPE::set_pr(double Pr1)
{
  Pr = Pr1;
}

int RTYPE::get_nadd()
{
  return nadd;
}

int RTYPE::get_nbrks()
{
  return nbrks;
}


void RTYPE::freemem()
{
  if (alloced!=1) return;

  delete [] add;
  delete [] brks;
  delete [] anames;
  delete [] anumbers;
  delete [] coordn;
  add = NULL;
  brks = NULL;
  anames = NULL;
  anumbers = NULL;
  coordn = NULL;

  natoms = 0;
  alloced = 0;

  return;
}

void RTYPE::init()
{
  inited = 1;
  alloced = 0;

  natoms = 0;
  add = NULL;
  brks = NULL;
  anames = NULL;
  anumbers = NULL;
  coordn = NULL;

  id = -1;

  return;
}
