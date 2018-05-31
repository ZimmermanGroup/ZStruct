// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "rxnftr.h"
using namespace std;


double RXNFTR::get_value_0(ICoord ic1, int nadd1, int* add1, int nbrk1, int* brk1)
{
 //counts how many atoms match input rxn (normalized by #atoms)

 //doesn't handle swap conditions

  double val = 0.;
//  if (nadd1+nbrk1<1) val = 0.;
  double val1;
  double val2;
 
 //at1/2 are attaching at1d/2d are attached

  at1.active = 1;
  at2.active = 1;
  for (int i=0;i<nadd1;i++)
  {
    val1 = val2 = 0.;

    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int e1 = ic1.anumbers[a1];
    int e2 = ic1.anumbers[a2];
    int c1 = ic1.coordn[a1];
    int c2 = ic1.coordn[a2];
    double q1 = ic1.q[a1];
    double q2 = ic1.q[a2];

    int nat1 = get_attached_all(ic1,a1,at1d);
    int nat2 = get_attached_all(ic1,a2,at2d);
    at1.set_element(e2); at1.set_coordn(c2); at1.set_q(q2);
    at2.set_element(e1); at2.set_coordn(c1); at2.set_q(q1);

    val1 = atom_match(e1,c1,q1,1,0,0,nat1,at1d,at1);
    val2 = atom_match(e2,c2,q2,1,0,0,nat2,at2d,at2);

    val += val1 + val2;
//    val += val1*val2;
  }

  at1.active = 0;
  at2.active = 0;
  for (int i=0;i<nbrk1;i++)
  {
    val1 = val2 = 0.;

    int b1 = brk1[2*i+0];
    int b2 = brk1[2*i+1];
    int e1 = ic1.anumbers[b1];
    int e2 = ic1.anumbers[b2];
    int c1 = ic1.coordn[b1];
    int c2 = ic1.coordn[b2];
    double q1 = ic1.q[b1];
    double q2 = ic1.q[b2];

    int nat1 = get_attached_all(ic1,b1,at1d);
    int nat2 = get_attached_all(ic1,b2,at2d);

    val1 = atom_match(e1,c1,q1,0,1,0,nat1,at1d,at1);
    val2 = atom_match(e2,c2,q2,0,1,0,nat2,at2d,at2);

    val += val1 + val2;
//    val += val1*val2;
  }


  for (int i=0;i<maxat;i++) 
  {
    at1d[i].active = 0;
    at2d[i].active = 0;
  }
  at1.active = 1;
  at2.active = 1;
  for (int i=0;i<nadd1;i++)
  for (int j=0;j<nbrk1;j++)
  {
    val1 = val2 = 0.;

    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int b1 = brk1[2*j+0];
    int b2 = brk1[2*j+1];

    int e1 = -1;
    int e2 = -1;
    int c1 = ic1.coordn[a1];
    int c2 = ic1.coordn[a2];
    int se1 = ic1.anumbers[a2]; //destination atoms
    int se2 = ic1.anumbers[a1];
    int sc1 = ic1.coordn[a2];
    int sc2 = ic1.coordn[a1];

    if (a1==b1)
      e1  = ic1.anumbers[a1];
    else if (a1==b2)
      e1  = ic1.anumbers[a1];
    if (a2==b2)
      e2  = ic1.anumbers[a2];
    else if (a2==b1)
      e2  = ic1.anumbers[a2];

    double q1 = ic1.q[a1];
    double q2 = ic1.q[a2];
    if (e1!=-1)
    {
      int nat1 = get_attached_all(ic1,a1,at1d);
      at1.set_element(se1); at1.set_coordn(sc1); at1.set_q(q2);
      val1 = atom_match(e1,c1,q1,0,0,1,nat1,at1d,at1);
    }
    if (e2!=-1)
    {
      int nat2 = get_attached_all(ic1,a2,at2d);
      at2.set_element(se2); at2.set_coordn(sc2); at2.set_q(q1);
      val2 = atom_match(e2,c2,q2,0,0,1,nat2,at2d,at2);
    }

    val += val1 + val2;
//    val += val1*val2;
  } //loop i,j over nadd1, nbrk1


  //printf(" rxnftrs value: %2.1f \n",val);

//  val = val / (nadd1 + nbrk1);

  //printf(" get_value: %4.1f \n",val);
  return val;
}


double RXNFTR::get_value_1(ICoord ic1, int nadd1, int* add1, int nbrk1, int* brk1)
{
 //doesn't handle swap conditions

  double val = 0.;
  if (nadd1+nbrk1<1) val = 0.;
  double val1;
  double val2;

  int nmatched = 0;
  int* matched = new int[natoms];
  for (int i=0;i<natoms;i++)
    matched[i] = 0;
 
#if 0
 //at1/2 are attaching at1d/2d are attached
  ATOM at1,at2,at1d,at2d;
  at1.reset();
  at2.reset();
  at1d.reset();
  at2d.reset();

  at1.active = 1;
  at2.active = 1;
  for (int i=0;i<nadd1;i++)
  {
    val1 = val2 = 0.;

    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    int e1 = ic1.anumbers[a1];
    int e2 = ic1.anumbers[a2];
    int c1 = ic1.coordn[a1];
    int c2 = ic1.coordn[a2];
    double q1 = ic1.q[a1];
    double q2 = ic1.q[a2];

  //handle just H for now
    if (e1==1) get_attached(ic1,a1,at1d); else at1d.active = 0;
    if (e2==1) get_attached(ic1,a2,at2d); else at1d.active = 0;
    at1.set_element(e2);
    at2.set_element(e1);

    val1 = atom_match(e1,c1,q1,1,0,0,nat1,at1d,at1,matched);
    val2 = atom_match(e2,c2,q2,1,0,0,nat2,at2d,at2,matched);

//    val += val1 + val2;
//    val += val1*val2;
  }

  at1.active = 0;
  at2.active = 0;
  for (int i=0;i<nbrk1;i++)
  {
    val1 = val2 = 0.;

    int b1 = brk1[2*i+0];
    int b2 = brk1[2*i+1];
    int e1 = ic1.anumbers[b1];
    int e2 = ic1.anumbers[b2];
    int c1 = ic1.coordn[b1];
    int c2 = ic1.coordn[b2];
    double q1 = ic1.q[b1];
    double q2 = ic1.q[b2];

  //handle just H for now
    if (e1==1) get_attached(ic1,b1,at1d); else at1d.active = 0;
    if (e2==1) get_attached(ic1,b2,at2d); else at2d.active = 0;

    val1 = atom_match(e1,c1,q1,0,1,0,nat1,at1d,at1,matched);
    val2 = atom_match(e2,c2,q2,0,1,0,nat2,at2d,at2,matched);

//    val += val1 + val2;
//    val += val1*val2;
  }

  //printf(" rxnftrs value: %2.1f \n",val);
#else
  printf(" WARNING: get_value_1 not implemented \n");
  exit(1);
#endif



 //return 1 if all atoms are matched
  for (int i=0;i<natoms;i++)
  if (matched[i])
    nmatched++;
  if (nmatched>=natoms)
    val = 1.;

  delete [] matched;

  //printf(" get_value: %4.1f \n",val);
  return val;
}


int RXNFTR::get_attached_all(ICoord ic1, int a1, ATOM* atd)
{
  //gets attached atoms, their elements + coordn
  //does not do abs

  //assumes atd's are init'd
  for (int i=0;i<maxat;i++)
    atd[i].active = 0;

  int nfound = 0;
  for (int i=0;i<ic1.nbonds;i++)
  {
    if (ic1.bonds[i][0]==a1)
    {
      int a2 = ic1.bonds[i][1];
      atd[nfound].set_element(ic1.anumbers[a2]);
      atd[nfound].set_coordn(ic1.coordn[a2]);
      atd[nfound++].active = 1;
    }
    else if (ic1.bonds[i][1]==a1)
    {
      int a2 = ic1.bonds[i][0];
      atd[nfound].set_element(ic1.anumbers[a2]);
      atd[nfound].set_coordn(ic1.coordn[a2]);
      atd[nfound++].active = 1;
    }
    if (nfound>=maxat)
      break;
  }
  //printf("   nfound: %i coordn: %i \n",nfound,ic1.coordn[a1]);

  return nfound;
}

void RXNFTR::get_attached(ICoord ic1, int a1, ATOM& atd1)
{
  atd1.active = 0;
  for (int i=0;i<ic1.nbonds;i++)
  {
    if (ic1.bonds[i][0]==a1)
    {
      int a2 = ic1.bonds[i][1];
      atd1.set_element(ic1.anumbers[a2]);
      atd1.set_coordn(ic1.coordn[a2]);
      atd1.active = 1;
      break;
    }
    else if (ic1.bonds[i][1]==a1)
    {
      int a2 = ic1.bonds[i][0];
      atd1.set_element(ic1.anumbers[a2]);
      atd1.set_coordn(ic1.coordn[a2]);
      atd1.active = 1;
      break;
    }
  }

  return;
}



int RXNFTR::count_parameters()
{
  int nparam = 0;
  for (int i=0;i<natoms;i++)
    nparam += atoms[i].count_parameters();

  p = nparam;
  return nparam;
}

int RXNFTR::get_element(int wa)
{
  return atoms[wa].get_element();
}

int RXNFTR::get_coordn(int wa)
{
  return atoms[wa].get_coordn();
}

double RXNFTR::get_q(int wa)
{
  return atoms[wa].get_q();
}

double RXNFTR::get_qsig(int wa)
{
  return atoms[wa].get_qsig();
}

void RXNFTR::get_abs(int wa, int& add1, int& brk1, int& swp1)
{
  atoms[wa].get_abs(add1,brk1,swp1);
  return;
}

void RXNFTR::reset(int natoms1)
{
  //reset atoms (blank)
  for (int i=0;i<natoms;i++)
    atoms[i].reset();
  for (int i=natoms;i<natoms1;i++)
    atoms[i].init();

  natoms = natoms1;
  new_atom = natoms - 1;

  return;
}

void RXNFTR::init(int natoms1) 
{
  max_atoms = 16;
  if (natoms1<1)
  {
    printf(" ERROR: cannot init with natoms < 1 \n");
    exit(-1);
  }
  if (natoms1>max_atoms)
  {
    printf(" ERROR: cannot init with natoms > %i \n",max_atoms);
  }
  quiet = 1;

  if (!quiet)
    printf("\n init'ing feature vector (%i) \n",natoms1);

  natoms = natoms1;
  atoms = new ATOM[max_atoms];
  new_atom = natoms1 - 1;

  //initialize atoms (blank)
  for (int i=0;i<natoms;i++)
    atoms[i].init();

 //for get_value
  at1.init(); at2.init();
  maxat = 8;
  at1d = new ATOM[maxat];
  at2d = new ATOM[maxat];
  for (int i=0;i<maxat;i++)
    at1d[i].init();
  for (int i=0;i<maxat;i++)
    at2d[i].init();


  return;
}

void RXNFTR::update(int natoms1)
{
  for (int i=natoms;i<natoms1;i++)
    atoms[i].init();
  for (int i=natoms;i>natoms1;i--)
    atoms[i-1].reset();
  natoms = natoms1;
  new_atom = natoms1 - 1;

  return;
}


void RXNFTR::freemem()
{
  for (int i=0;i<natoms;i++)
    atoms[i].freemem();
  delete [] atoms;

  return;
}

void RXNFTR::print_atoms()
{
  printf("   atoms: \n");
  for (int i=0;i<natoms;i++)
    atoms[i].print();

  return;
}

void RXNFTR::set_element(int wa, int e1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_element(e1);
  return;
}

void RXNFTR::set_coordn(int wa, int c1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_coordn(c1);
  return;
}

void RXNFTR::set_q(int wa, double q1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_q(q1);
  return;
}

void RXNFTR::set_qsig(int wa, double qsig1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_qsig(qsig1);
  return;
}

void RXNFTR::set_abs(int wa, int add1, int brk1, int swp1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_abs(add1,brk1,swp1);
  return;
}

void RXNFTR::set_attached(int wa, ATOM at1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attached(at1);
}

void RXNFTR::set_attached(int wa, int e1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attached(e1);
}

void RXNFTR::set_attached(int wa, int e1, int c1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attached(e1,c1);
}

void RXNFTR::set_attached(int wa, int e1, int c1, double q1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attached(e1,c1,q1);
}

void RXNFTR::set_attaching(int wa, ATOM at2)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attaching(at2);
}

void RXNFTR::set_attaching(int wa, int e1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attaching(e1);
}

void RXNFTR::set_attaching(int wa, int e1, int c1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attaching(e1,c1);
}

void RXNFTR::set_attaching(int wa, int e1, int c1, double q1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].set_attaching(e1,c1,q1);
}

void RXNFTR::delete_attached(int wa)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].delete_attached();
}

void RXNFTR::delete_attaching(int wa)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].delete_attaching();
}

int RXNFTR::get_nattached(int wa)
{
  if (wa>=natoms) return 0;
  if (wa==-1) wa = new_atom;
  return atoms[wa].get_nattached();
}

int RXNFTR::get_nattaching(int wa)
{
  if (wa>=natoms) return 0;
  if (wa==-1) wa = new_atom;
  return atoms[wa].get_nattaching();
}

void RXNFTR::get_attached(int wa, ATOM& at1)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].get_attached(at1);
}

void RXNFTR::get_attaching(int wa, ATOM& at2)
{
  if (wa>=natoms) return;
  if (wa==-1) wa = new_atom;
  atoms[wa].get_attaching(at2);
}

double RXNFTR::atom_match(int e1, int c1, double q1, int add1, int brk1, int swp1, int nat1, ATOM* at1, ATOM at2)
{
  double val = 0;

  for (int i=0;i<natoms;i++)
  {
    val += atoms[i].compare(e1,c1,q1,add1,brk1,swp1,nat1,at1,at2); 
  }
  val = val / natoms;

  //printf(" atom_match val: %2.1f \n",val);
  return val;
}

double RXNFTR::atom_match(int e1, int c1, double q1, int add1, int brk1, int swp1, int nat1, ATOM* at1, ATOM at2, int* matched)
{
  double val = 0;

 //label all matching atoms
  for (int i=0;i<natoms;i++)
  {
    val = atoms[i].compare(e1,c1,q1,add1,brk1,swp1,nat1,at1,at2); 
  //CPMZ may need to change this criterion
    if (val>0.) matched[i] = 1;
  }

  //printf(" atom_match val: %2.1f \n",val);
  return val;
}

