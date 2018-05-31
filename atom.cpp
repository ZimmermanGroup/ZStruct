// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "atom.h"
using namespace std;

//NOTE: remove init'd error messages

//attached and attaching match are binary
//compare now returns double



int ATOM::compare(int e1, int c1, int add1, int brk1, int swp1)
{
  int same = 1;

  if      (!element_match(e1)) same = 0;
  else if (!coordn_match(c1)) same = 0;
  else if (!add_match(add1)) same = 0;
  else if (!brk_match(brk1)) same = 0;
  else if (!swp_match(swp1)) same = 0;

  return same;
}

double ATOM::compare(int e1, int c1, double q1, int add1, int brk1, int swp1)
{
  int same = 1;

  if      (!element_match(e1)) same = 0;
  else if (!coordn_match(c1)) same = 0;
  else if (!add_match(add1)) same = 0;
  else if (!brk_match(brk1)) same = 0;
  else if (!swp_match(swp1)) same = 0;

  double dsame = same;
  if (same)
    dsame = q_eval(q1);

  return dsame;
}

double ATOM::compare(int e1, int c1, double q1, int add1, int brk1, int swp1, int nat1, ATOM* at1, ATOM at2)
{
  int same = compare(e1,c1,add1,brk1,swp1);

//  return same;

  if (same && nat1>0)
  {
    int found = 0;
    for (int i=0;i<nat1;i++)
    if (attached_match(at1[i]))
      found = 1;
    if (!found)  same = 0;
  }

 //q evaluated only for attaching
  double dsame = same;
  if (same)
  {
    dsame = attaching_match(at2);
    dsame *= q_eval(q1);
  }

#if 0
  if (qsig<1000. && same)
    printf(" q1: %4.2f q: %4.2f dsame: %4.3f \n",q1,q,dsame);
#endif

  return dsame;
}

double ATOM::q_eval(double q1)
{
  if (qsig>=1000. || q1==0.0) return 1.;

  double qe;
  double d = q1 - q;
//  qe = exp(-d*d/qsig);
  qe = exp(-fabs(d)/qsig);

  return qe;
}

void ATOM::print()
{
  printf("    element: %2i coordn: %2i q: %4.1f ABS: %2i %2i %2i",elem,coordn,q,add,brk,swp);
  if (qsig<1000.) printf(" *");
  else printf("  ");

  for (int i=0;i<nattached;i++)
  {
    printf(" attached to:  %i (c%i",attached[i].get_element(),attached[i].get_coordn());
    if (attached[i].get_qsig()<1000.)
      printf(" %4.1f",attached[i].get_q());
    printf(")");
  }
  for (int i=0;i<nattaching;i++)
  {
    printf(" attaching to: %i (c%i",attaching[i].get_element(),attaching[i].get_coordn());
    if (attaching[i].get_qsig()<1000.)
      printf(" %4.1f",attaching[i].get_q());
    printf(")");
  }
#if 0
  for (int i=0;i<nattached;i++)
    attached[i].print();
  for (int i=0;i<nattaching;i++)
    attaching[i].print();
#endif

  printf("\n");
}

void ATOM::freemem()
{
#if 0
  if (inited!=1)
  {
    printf(" ERROR in ATOM::freemem, not init'd \n");
    exit(1);
  }
#endif

  if (attached!=NULL)
  {
    for (int i=0;i<nattached;i++)
      attached[i].freemem();
    delete [] attached;
    attached = NULL;
  }
  if (attaching!=NULL)
  {
    for (int i=0;i<nattaching;i++)
      attaching[i].freemem();
    delete [] attaching;
    attaching = NULL;
  }
  nattached = 0;
  nattaching = 0;

  inited = 0;
}

void ATOM::reset()
{
#if 0
  if (inited!=1)
  {
    printf(" ERROR in ATOM::reset, not init'd \n");
    exit(1);
  }
#endif

  add = -1;
  brk = -1;
  swp = -1;

  elem = -1;
  coordn = -1;

  if (attached!=NULL)
  {
    for (int i=0;i<nattached;i++)
      attached[i].freemem();
    delete [] attached;
    attached = NULL;
  }
  if (attaching!=NULL)
  {
    for (int i=0;i<nattaching;i++)
      attaching[i].freemem();
    delete [] attaching;
    attaching = NULL;
  }

  nattached = 0;
  nattaching = 0;

  return;
} 

void ATOM::set_element(int e1)
{
  elem = e1;
  return;
}

void ATOM::set_coordn(int c1)
{
  coordn = c1;
  return;
}

void ATOM::set_q(double q1)
{
  q = q1;
  return;
}

void ATOM::set_qsig(double qsig1)
{
  qsig = qsig1;
  return;
}

void ATOM::set_abs(int add1, int brk1, int swp1)
{
  add = add1;
  brk = brk1;
  swp = swp1;
  return;
}

void ATOM::set_attached(int e1)
{
 //only implemented for one attached atom
#if 0
  if (inited!=1)
  {
    printf(" ERROR in ATOM::freemem, not init'd \n");
    exit(1);
  }
#endif

  if (attached!=NULL)
  {
    for (int i=0;i<nattached;i++)
      attached[i].freemem();
    delete [] attached; 
    attached = NULL;
  }

  nattached = 1;
  if (attached==NULL)
  {
    attached = new ATOM[nattached];
    for (int i=0;i<nattached;i++)
      attached[i].init();
  }

  attached[0].set_element(e1);

  //not yet using coordn or abs
//  attached[0].set_coordn(at1.get_coordn());
//  attached[0].set_abs(-1,-1,-1);

  return;
}

void ATOM::set_attached(int e1, int c1)
{
  set_attached(e1);
  attached[0].set_coordn(c1);

  return;
}

void ATOM::set_attached(int e1, int c1, double q1)
{
  set_attached(e1);
  attached[0].set_coordn(c1);
  attached[0].set_q(q1);
  attached[0].set_qsig(0.4);

  return;
}


void ATOM::set_attached(ATOM at1)
{
 //only implemented for one attached atom

  if (attached!=NULL)
  {
    for (int i=0;i<nattached;i++)
      attached[i].freemem();
    delete [] attached;
    attached = NULL;
  }

  nattached = 1;
  if (attached==NULL)
  {
    attached = new ATOM[nattached];
    for (int i=0;i<nattached;i++)
      attached[i].init();
  }

  attached[0].set_element(at1.get_element());

  //not yet using coordn or abs
//  attached[0].set_coordn(at1.get_coordn());
//  attached[0].set_abs(-1,-1,-1);


  return;
}

void ATOM::set_attaching(int e1)
{
 //only implemented for one attached atom

  if (attaching!=NULL)
  {
    for (int i=0;i<nattaching;i++)
      attaching[i].freemem();
    delete [] attaching;
    attaching = NULL;
  }

  nattaching = 1;
  if (attaching==NULL)
  {
    attaching = new ATOM[nattaching];
    for (int i=0;i<nattaching;i++)
      attaching[i].init();
  }

  attaching[0].set_element(e1);

  //not yet using coordn or abs
//  attaching[0].set_coordn(at2.get_coordn());
//  attaching[0].set_abs(-1,-1,-1);


  return;
}

void ATOM::set_attaching(int e1, int c1)
{
  set_attaching(e1);
  attaching[0].set_coordn(c1);

  return;
}

void ATOM::set_attaching(int e1, int c1, double q1)
{
  set_attaching(e1);
  attaching[0].set_coordn(c1);
  attaching[0].set_q(q1);
  attaching[0].set_qsig(0.4);

  return;
}

void ATOM::set_attaching(ATOM at2)
{
 //only implemented for one attached atom

  if (attaching!=NULL)
  {
    for (int i=0;i<nattaching;i++)
      attaching[i].freemem();
    delete [] attaching;
    attaching = NULL;
  }

  nattaching = 1;
  if (attaching==NULL)
  {
    attaching = new ATOM[nattaching];
    for (int i=0;i<nattaching;i++)
      attaching[i].init();
  }

  attaching[0].set_element(at2.get_element());

  //not yet using coordn or abs
//  attaching[0].set_coordn(at2.get_coordn());
//  attaching[0].set_abs(-1,-1,-1);


  return;
}

void ATOM::delete_attached()
{
  if (attached!=NULL)
  {
    for (int i=0;i<nattached;i++)
      attached[i].freemem();
    delete [] attached;
    attached = NULL;
  }

  nattached = 0;

  return;
}

void ATOM::delete_attaching()
{
  if (attaching!=NULL)
  {
    for (int i=0;i<nattaching;i++)
      attaching[i].freemem();
    delete [] attaching;
    attaching = NULL;
  }

  nattaching = 0;

  return;
}


int ATOM::get_element()
{
  return elem;
}

int ATOM::get_coordn()
{
  return coordn;
}

double ATOM::get_q()
{
  return q;
}

double ATOM::get_qsig()
{
  return qsig;
}

void ATOM::get_abs(int& add1, int& brk1, int& swp1) 
{
  add1 = add;
  brk1 = brk;
  swp1 = swp;
  return;
}

int ATOM::get_nattaching()
{
  return nattaching;
}

int ATOM::get_nattached()
{
  return nattached;
}

void ATOM::get_attached(ATOM& at1)
{
  if (nattached<1) return;

  int add1,brk1,swp1;
  int e1 = attached[0].get_element();
  int c1 = attached[0].get_coordn();
  attached[0].get_abs(add1,brk1,swp1);

  at1.set_element(e1);
  at1.set_coordn(c1);
  at1.set_abs(add1,brk1,swp1);

  return;
}

void ATOM::get_attaching(ATOM& at1)
{
  if (nattaching<1) return;

  int add1,brk1,swp1;
  int e1 = attaching[0].get_element();
  int c1 = attaching[0].get_coordn();
  attaching[0].get_abs(add1,brk1,swp1);

  at1.set_element(e1);
  at1.set_coordn(c1);
  at1.set_abs(add1,brk1,swp1);

  return;
}


int ATOM::element_match(int e1)
{
  if (elem==-1 || e1==-1) return 1;

  int found = 0;
  if (e1==elem)
    found = 1;

  return found;
}

int ATOM::coordn_match(int c1)
{
  if (coordn==-1 || c1==-1) return 1;

  int found = 0;
  if (c1==coordn)
    found = 1;

  return found;
}

int ATOM::add_match(int add1)
{
  if (add==-1 || add1==-1) return 1;

  int found = 0;
  if (add1==add)
    found = 1;

  return found;
}

int ATOM::brk_match(int brk1)
{
  if (brk==-1 || brk1==-1) return 1;

  int found = 0;
  if (brk1==brk)
    found = 1;

  return found;
}

int ATOM::swp_match(int swp1)
{
  if (swp==-1 || swp1==-1) return 1;

  int found = 0;
  if (swp1==swp)
    found = 1;

  return found;
}


int ATOM::attached_match(ATOM at1)
{
  if (nattached<1 || !at1.active)
    return 1;

  int found = 0;
  int e1,c1,add1,brk1,swp1;
  e1 = at1.get_element();
  c1 = at1.get_coordn();
  //double q1 = at1.get_q();
  at1.get_abs(add1,brk1,swp1);
  for (int i=0;i<nattached;i++)
  {
    found += attached[i].compare(e1,c1,add1,brk1,swp1);
  }
  //printf(" ad%i",found); 

  return found;
}

double ATOM::attaching_match(ATOM at2)
{
  if (nattaching<1 || !at2.active)
    return 1.;

  double found = 0;
  int e1,c1,add1,brk1,swp1;
  e1 = at2.get_element();
  c1 = at2.get_coordn();
  double q1 = at2.get_q();
  at2.get_abs(add1,brk1,swp1);
  for (int i=0;i<nattaching;i++)
  {
    found += attaching[i].compare(e1,c1,q1,add1,brk1,swp1);
  }

  //at2.print();
//  if (attaching[0].get_qsig()<1000. && found>0.01)
//    printf(" q1: %4.1f q: %4.1f found: %4.2f \n",q1,attaching[0].get_q(),found);
//    printf(" e1: %2i c1: %2i q1: %4.1f q: %4.1f found: %4.2f \n",e1,c1,q1,attaching[0].get_q(),found);

  return found;
}


int ATOM::count_parameters()
{
  int nparam = 0;

  if (elem>-1) nparam++;
  if (coordn>-1) nparam++;
  if (qsig<1000.) nparam++;
  if (add>-1) nparam++;
  if (brk>-1) nparam++;
  if (swp>-1) nparam++;

#if 1
  for (int i=0;i<nattached;i++)
    nparam += attached[i].count_parameters();
  for (int i=0;i<nattaching;i++)
    nparam += attaching[i].count_parameters();
#endif

  p = nparam;
  return nparam;
}

void ATOM::init()
{
  active = 0;

  add = -1;
  brk = -1;
  swp = -1;

  elem = -1;
  coordn = -1;

  q = 0.;
  qsig = 1000.;

  nattached = 0;
  nattaching = 0;

 //for later
  attached = NULL;
  attaching = NULL;

  inited = 1;

  return;
}
