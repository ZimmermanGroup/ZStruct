// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "rxndb.h"
using namespace std;

///NOTE: now 60kT 
#define USE_KNNR 1
#define KNNR_K 4
#define DO_NOT_ADD_DUPLICATES 1

#define TEST_GA 1
#define MAX_N 1000
//for KNNR -- 686 for 1b? 

#define TRAINSHIFT 15.0
#define PTHRESH0 0.4
#define LAMBDA 0.00
//no evidence NULL_ATOMS is helpful
#define NULL_ATOMS 0

//GA: need to implement pthresh calc
//GA: replace bad features with A-B connections?
//NOTE: charges on attaching atoms OFF


void RXNDB::make_decision_tree()
{
  if (ga_inited<1) init_ga(1);

  printf("\n creating simple decision tree \n\n");


  naddbrktypes = 12; //a1, b1, a2, b2, a1b1, a2b1, a1b2, a2b2, a3b2, a2b3, a3b3, axbx
  if (addbrk_layer!=NULL) delete [] addbrk_layer;
  addbrk_layer = new double[naddbrktypes];
  for (int i=0;i<naddbrktypes;i++) addbrk_layer[i] = 0.;

  if (allow_add_brk_elem!=NULL) delete [] allow_add_brk_elem;
  allow_add_brk_elem = new double[nelemf*nelemf];
  for (int i=0;i<nelemf*nelemf;i++) allow_add_brk_elem[i] = 0.;

  int maxc = 1;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int c = climit_h[e1];
    if (c>maxc)
      maxc = c;
  }
  maxc++;

  double* weightc = new double[nelemf*maxc];
  double* weightp = new double[nelemf*nelemf];
  for (int i=0;i<nelemf*maxc;i++) weightc[i] = 0.;
  for (int i=0;i<nelemf*nelemf;i++) weightp[i] = 0.;

  double* weightaddc = new double[nelemf*maxc];
  double* weightpadd = new double[nelemf*nelemf];
  double* weightbrkc = new double[nelemf*maxc];
  double* weightpbrk = new double[nelemf*nelemf];
  for (int i=0;i<nelemf*maxc;i++) weightaddc[i] = 0.;
  for (int i=0;i<nelemf*nelemf;i++) weightpadd[i] = 0.;
  for (int i=0;i<nelemf*maxc;i++) weightbrkc[i] = 0.;
  for (int i=0;i<nelemf*nelemf;i++) weightpbrk[i] = 0.;

  int naddmax = 0;
  int nbrkmax = 0;
  for (int i=0;i<N;i++)
  {
    int wr = i;
    int natoms1 = natoms[wr];
    string* anames1 = anames[wr];
    int* anumbers1 = anumbers[wr];
    int* coordnr1 = coordnr[wr];
    double* xyzr1 = xyzr[wr];
    int nadd1 = nadd[wr];
    int* add1 = add[wr];
    int nbrks1 = nbrks[wr];
    int* brks1 = brks[wr];

    if (naddmax<nadd1) naddmax = nadd1;
    if (nbrkmax<nbrks1) nbrkmax = nbrks1;


    double val = Pr[i];
    int abtype = get_abtype(i,nadd1,nbrks1);

   //nadd/nbrk possible
    addbrk_layer[abtype] += val;

    if (abtype==11 && val>0.5)  printf("  big mover (add: %i brk: %i): %4i/%4i \n",nadd1,nbrks1,i,ids[i]);
    //printf(" i: %i Pr: %4.3f \n",i,val);

   //add/brk possible at element combinations
    for (int j=0;j<nadd1;j++)
    {
      int a1 = add1[2*j+0];
      int a2 = add1[2*j+1];
      int e1 = anumbers1[a1];
      int e2 = anumbers1[a2];
      int er1 = remap[e1];
      int er2 = remap[e2];
      int c1 = coordnr1[a1];
      int c2 = coordnr1[a2];

      weightc[er1*maxc+c1] += val; weightc[er2*maxc+c2] += val;
      weightp[er1*nelemf+er2] += val;
      if (er1!=er2)
        weightp[er2*nelemf+er1] += val;

      weightaddc[er1*maxc+c1] += val; weightaddc[er2*maxc+c2] += val;
      weightpadd[er1*nelemf+er2] += val;
      if (er1!=er2)
        weightpadd[er2*nelemf+er1] += val;
    }
    for (int j=0;j<nbrks1;j++)
    {
      int b1 = brks1[2*j+0];
      int b2 = brks1[2*j+1];
      int e1 = anumbers1[b1];
      int e2 = anumbers1[b2];
      int er1 = remap[e1];
      int er2 = remap[e2];
      int c1 = coordnr1[b1];
      int c2 = coordnr1[b2];

      weightc[er1*maxc+c1] += val; weightc[er2*maxc+c2] += val;
      weightp[er1*nelemf+er2] += val;
      if (er1!=er2)
        weightp[er2*nelemf+er1] += val;

      weightbrkc[er1*maxc+c1] += val; weightbrkc[er2*maxc+c2] += val;
      weightpbrk[er1*nelemf+er2] += val;
      if (er1!=er2)
        weightpbrk[er2*nelemf+er1] += val;
    }
  } //loop i over N data points



  
  for (int i=0;i<nelemf*maxc;i++) weightc[i] = weightc[i] / N;
  for (int i=0;i<nelemf*nelemf;i++) weightp[i] = weightp[i] / N;
  for (int i=0;i<nelemf*maxc;i++) weightaddc[i] = weightaddc[i] / N;
  for (int i=0;i<nelemf*nelemf;i++) weightpadd[i] = weightpadd[i] / N;
  for (int i=0;i<nelemf*maxc;i++) weightbrkc[i] = weightbrkc[i] / N;
  for (int i=0;i<nelemf*nelemf;i++) weightpbrk[i] = weightpbrk[i] / N;

  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    string anam = PTable::atom_name(e1);

    for (int j=climit_l[e1];j<climit_h[e1]+1;j++)
      printf(" element %s with coordn %i has weights: %4.3f  %4.3f  %4.3f \n",anam.c_str(),j,weightc[i*maxc+j],weightaddc[i*maxc+j],weightbrkc[i*maxc+j]);
  }

  printf("\n pair weights: \n ");
  for (int i=0;i<nelemf;i++)
    printf(" %5s",PTable::atom_name(emap[i]).c_str());
  printf("\n");
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    string anam1 = PTable::atom_name(e1);
    printf(" %2s:",anam1.c_str());
    for (int j=0;j<=i;j++)
    {
      //int e2 = emap[j];
      //string anam2 = PTable::atom_name(e2);

      printf(" %4.3f",weightp[i*nelemf+j]);
    }
    printf("\n");
  }
  printf("\n add weights: \n ");
  for (int i=0;i<nelemf;i++)
    printf(" %5s",PTable::atom_name(emap[i]).c_str());
  printf("\n");
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    string anam1 = PTable::atom_name(e1);
    printf(" %2s:",anam1.c_str());
    for (int j=0;j<=i;j++)
    {
      //int e2 = emap[j];
      //string anam2 = PTable::atom_name(e2);

      printf(" %4.3f",weightpadd[i*nelemf+j]);
    }
    printf("\n");
  }
  printf("\n break weights: \n ");
  for (int i=0;i<nelemf;i++)
    printf(" %5s",PTable::atom_name(emap[i]).c_str());
  printf("\n");
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    string anam1 = PTable::atom_name(e1);
    printf(" %2s:",anam1.c_str());
    for (int j=0;j<=i;j++)
    {
      //int e2 = emap[j];
      //string anam2 = PTable::atom_name(e2);

      printf(" %4.3f",weightpbrk[i*nelemf+j]);
    }
    printf("\n");
  }



  for (int i=0;i<nelemf;i++)
  for (int j=0;j<nelemf;j++)
  if (weightp[i*nelemf+j]>0.001) //CPMZ parameter
    allow_add_brk_elem[i*nelemf+j] = 1.0;

  printf("\n allow_add_brk_elem: \n    ");
  for (int i=0;i<nelemf;i++)
    printf(" %s",PTable::atom_name(emap[i]).c_str());
  printf("\n");
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    string anam1 = PTable::atom_name(e1);
    printf(" %2s:",anam1.c_str());
    for (int j=0;j<=i;j++)
      printf(" %1.0f",allow_add_brk_elem[i*nelemf+j]);
    printf("\n");
  }

  printf("\n                add1    brk1    add2    brk2    a1b1    a2b1    a1b2    a2b2    a3b2    a2b3    a3a3    axbx \n"); 
  printf(" addbrk_layer:");
  for (int i=0;i<naddbrktypes;i++)
    printf(" %7.3f",addbrk_layer[i]);
  printf("\n");

  delete [] weightc;
  delete [] weightp;
  delete [] weightaddc;
  delete [] weightpadd;
  delete [] weightbrkc;
  delete [] weightpbrk;


  make_combo_tree(naddmax,nbrkmax);



  printf("\n\n\n testing the Trees \n");
  //test trees
  int npass = 0;
  for (int i=0;i<N;i++)   
  {
    int wr = i;
    int natoms1 = natoms[wr];
    string* anames1 = anames[wr];
    int* anumbers1 = anumbers[wr];
    int* coordnr1 = coordnr[wr];
    double* xyzr1 = xyzr[wr];
    int nadd1 = nadd[wr];
    int* add1 = add[wr];
    int nbrks1 = nbrks[wr];
    int* brks1 = brks[wr];
    double* qrc1 = qrc[wr];

    double val = Pr[wr];
    int abtype = get_abtype(wr,nadd1,nbrks1);

    int tree0 = tree_screen_0(nadd1,nbrks1);
    double tree1 = tree_screen_1(natoms1,anames1,anumbers1,xyzr1,nadd1,add1,nbrks1,brks1,qrc1);
    double tree2 = tree_screen_2(natoms1,anames1,anumbers1,coordnr1,xyzr1,nadd1,add1,nbrks1,brks1,qrc1);

    printf(" rxn %i (ids: %4i nadd: %i nbrks: %i) tree0: %i tree1: %4.3f tree2: %4.3f  \n",wr,ids[wr],nadd1,nbrks1,tree0,tree1,tree2);

    if (tree0 && tree1 > 0.01 && tree2 > 0.01) npass++;

  }

  printf("\n pass: %i of %i fraction: %4.3f \n",npass,N,1.*npass/N);

printf("\n ending early ! \n");
exit(1);

  return;
}


void RXNDB::make_combo_tree(int naddmax, int nbrkmax)
{
//not using naddmax, nbrkmax?

  printf("\n\n in make_combo_tree \n");


  if (trxns!=NULL) trxns = NULL;
  trxns = new RTYPE*[naddbrktypes];

  if (nrtypes!=NULL) delete [] nrtypes;
  nrtypes = new int[naddbrktypes];
  for (int i=0;i<naddbrktypes;i++) nrtypes[i] = 0;

  int* maxabtypes = new int[naddbrktypes];
  for (int i=0;i<naddbrktypes;i++) maxabtypes[i] = 0;


 //count max # of each type for allocation
  for (int i=0;i<N;i++)
  {
    int wr = i;
    int nadd1 = nadd[wr];
    int nbrks1 = nbrks[wr];

    int abtype = get_abtype(i,nadd1,nbrks1);

    maxabtypes[abtype]++;
  }

  for (int i=0;i<naddbrktypes;i++)
    trxns[i] = new RTYPE[maxabtypes[i]];


  int maxchgatoms = 20;
  int* anumbers1 = new int[maxchgatoms];
  string* anames1 = new string[maxchgatoms];
  int* coordn1 = new int[maxchgatoms];
  int* add1 = new int[maxchgatoms];
  int* brks1 = new int[maxchgatoms];
  for (int i=0;i<maxchgatoms;i++)
    add1[i] = brks1[i] = i;

 //get all unique of each naddbrktype
  for (int i=0;i<N;i++)
  {
    int wr = i;
    string* anames0 = anames[wr];
    int* anumbers0 = anumbers[wr];
    int* coordn0 = coordnr[wr];
    int nadd1 = nadd[wr];
    int* add0 = add[wr];
    int nbrks1 = nbrks[wr];
    int* brks0 = brks[wr];
    int natoms1 = 2*nadd1 + 2*nbrks1;
    if (natoms1>maxchgatoms) 
    {
      printf(" ERROR: too many changes in rxn %i/%i \n",i,ids[i]);
      exit(1);
    }

    int naf = 0;
    for (int j=0;j<2*nadd1;j++)
    {
      int a1 = add0[j];
      coordn1[naf] = coordn0[a1]; //alt: set to -1
      coordn1[naf] = -1;
      anames1[naf] = anames0[a1]; //alt: set to X
//      anames1[naf] = "X";
      anumbers1[naf] = anumbers0[a1]; //alt: set to -1
//      anumbers1[naf] = -1;
     //add1[j] = naf; //set by default
      naf++;
    }
    for (int j=0;j<2*nbrks1;j++)
    {
      int b1 = brks0[j];
      coordn1[naf] = coordn0[b1];
      coordn1[naf] = -1;
      anames1[naf] = anames0[b1];
//      anames1[naf] = "X";
      anumbers1[naf] = anumbers0[b1];
//      anumbers1[naf] = -1;
      brks1[j] = naf;
      naf++;
    }

    if (naf!=natoms1)
    {
      printf("\n naf not equal to natoms1: %i %i \n",naf,natoms1);
      exit(1);
    }


    double val = Pr[i];
    int abtype = get_abtype(i,nadd1,nbrks1);
    int nrt = nrtypes[abtype];

    int found = -1;
    for (int j=0;j<nrt;j++)
    if (trxns[abtype][j].match(natoms1,anumbers1,coordn1,nadd1,add1,nbrks1,brks1))
    {
      //printf("   match found: %i \n",j);
      found = j;
      break;
    }
    if (found==-1)
    {
      //printf("  new rxn (%4i): natoms %i nadd %i nbrks %i \n",ids[i],natoms1,nadd1,nbrks1);
      trxns[abtype][nrt].init();
      trxns[abtype][nrt].id = ids[i];
      trxns[abtype][nrt].set_pr(val); 
      trxns[abtype][nrt].set_rxn(natoms1,anames1,anumbers1,coordn1,nadd1,add1,nbrks1,brks1);
      nrtypes[abtype]++;
    }
    else if (val>trxns[abtype][found].get_pr())
    {
      //printf("   is higher Pr than previous %4.3f/%4.3f \n",val,trxns[abtype][found].get_pr());
      trxns[abtype][nrt].id = ids[i];
      trxns[abtype][nrt].set_pr(val);
    }
    //else
    //  printf("  aleady found higher Pr than %4.3f/%4.3f \n",val,trxns[abtype][found].get_pr());

    

  } //loop i over N


  printf("\n printing rtypes: \n");
  for (int i=0;i<naddbrktypes;i++)
  if (nrtypes[i]>0)
  {
    int nadd0 = trxns[i][0].get_nadd();
    int nbrks0 = trxns[i][0].get_nbrks();
    printf("\n  type %i is nadd %i nbrk %i \n",i,nadd0,nbrks0);

    for (int j=0;j<nrtypes[i];j++)
    {
//      double val = trxns[i][j].get_pr();
//      printf("    Pr: %4.3f \n",val);
      trxns[i][j].print();
    }
  }
  else
    printf("\n  no data for type %i \n",i);


  delete [] add1;
  delete [] brks1;
  delete [] anames1;
  delete [] anumbers1;
  delete [] coordn1;
  delete [] maxabtypes;

#if 0
  printf("\n\n Now printing entire data set \n");
  print_reactions_2();
#endif

  return;
}




int RXNDB::tree_screen_0(int nadd0, int nbrks0)
{
  if (addbrk_layer==NULL || naddbrktypes<1) return 1;

  double val;
  int abtype = get_abtype(0,nadd0,nbrks0);
  val = addbrk_layer[abtype];

  int pass = 0;
  if (val>0.001)
    pass = 1;

  return pass;
}


double RXNDB::tree_screen_1(int natoms0, string* anames0, int* anumbers0, double* xyzr0, int nadd0, int* add0, int nbrks0, int* brks0, double* qrc0)
{
  double val = 1.;

  if (allow_add_brk_elem==NULL)
    return 1.0;

  //take product of allowed pairings..
  for (int i=0;i<nadd0;i++)
  {
    int a1 = add0[2*i+0];
    int a2 = add0[2*i+1];
    int e1 = anumbers0[a1];
    int e2 = anumbers0[a2];
    int e1m = remap[e1];
    int e2m = remap[e2];
    string anam1 = PTable::atom_name(e1);
    string anam2 = PTable::atom_name(e2);

    val *= allow_add_brk_elem[e1m*nelemf+e2m];
    //printf("  adding %s %s allowed: %2.1f \n",anam1.c_str(),anam2.c_str(),allow_add_brk_elem[e1m*nelemf+e2m]);
  }

  for (int i=0;i<nbrks0;i++)
  {
    int b1 = brks0[2*i+0];
    int b2 = brks0[2*i+1];
    int e1 = anumbers0[b1];
    int e2 = anumbers0[b2];
    int e1m = remap[e1];
    int e2m = remap[e2];
    string anam1 = PTable::atom_name(e1);
    string anam2 = PTable::atom_name(e2);

    val *= allow_add_brk_elem[e1m*nelemf+e2m];
    //printf("  brking %s %s allowed: %2.1f \n",anam1.c_str(),anam2.c_str(),allow_add_brk_elem[e1m*nelemf+e2m]);
  }

  return val;
}


double RXNDB::tree_screen_2(int natoms0, string* anames0, int* anumbers0, int* coordn0, double* xyzr0, int nadd0, int* add0, int nbrks0, int* brks0, double* qrc0)
{
  double val = 1.;

  int maxchgatoms = 20;
  int nadd1 = nadd0;
  int nbrks1 = nbrks0;
  int* anumbers1 = new int[maxchgatoms];
  string* anames1 = new string[maxchgatoms];
  int* coordn1 = new int[maxchgatoms];
  int* add1 = new int[maxchgatoms];
  int* brks1 = new int[maxchgatoms];
  for (int i=0;i<maxchgatoms;i++)
    add1[i] = brks1[i] = i;

  int naf = 0;
  for (int j=0;j<2*nadd0;j++)
  {
    int a1 = add0[j];
    coordn1[naf] = coordn0[a1]; 
    anames1[naf] = anames0[a1]; 
    anumbers1[naf] = anumbers0[a1]; 
   //add1[j] = naf; //set by default
    naf++;
  }
  for (int j=0;j<2*nbrks0;j++)
  {
    int b1 = brks0[j];
    coordn1[naf] = coordn0[b1];
    anames1[naf] = anames0[b1];
    anumbers1[naf] = anumbers0[b1];
    brks1[j] = naf;
    naf++;
  }

  int natoms1 = naf;
  int abtype = get_abtype(nadd0,nbrks0);
  int nrt = nrtypes[abtype];

  int found = -1;
  for (int j=0;j<nrt;j++)
  if (trxns[abtype][j].match(natoms1,anumbers1,coordn1,nadd1,add1,nbrks1,brks1))
  {
    //printf("   match found: %i \n",j);
    found = j;
    break;
  }
  if (found==-1)
  {
    printf("  no info on this reaction! \n");
  }
  else
  {
    printf("  found match \n");
    trxns[abtype][found].print();

    val = trxns[abtype][found].get_pr();
  }


  delete [] add1;
  delete [] brks1;
  delete [] anumbers1;
  delete [] anames1;
  delete [] coordn1;

  return val;
}



void RXNDB::element_analysis(int wg)
{
  if (wg<0)
    wg = best_fv[0];

  printf("\n\n now analyzing elemental reactivity (wg: %i) \n",wg);
  if (p[wg]<1) 
  {
    printf("  cannot make estimate, p=0 \n");
    return;
  }

#if 1
  ICoord ic1;
  ic1.alloc(3);

  int* anumbers1 = new int[3];
  string* anames1 = new string[3];
  double* xyz1 = new double[9];
  for (int i=0;i<9;i++) xyz1[i] = 0.;
  xyz1[3] =  3.0;
  xyz1[6] = -1.0;
  double* qrc1 = new double[3];
  for (int i=0;i<3;i++) qrc1[i] = 0.;
  ic1.q = qrc1;

  int* add1 = new int[4];
  add1[0] = 0; add1[1] = 1;
  int* brk1 = new int[4];
  brk1[0] = 0; brk1[1] = 1;

  double* X1 = new double[p[wg]+1];
  double* X2 = new double[p[wg]+1];
  double* X3 = new double[p[wg]+1];
  for (int i=0;i<nelemf;i++)
//  for (int j=0;j<=i;j++)
  {
    int e1 = anumbers1[0] = emap[i]; //reacting atom
//    int e2 = anumbers1[1] = emap[j]; //adding to
    int e2 = anumbers1[1] = 1;
    anames1[0] = PTable::atom_name(e1);
//    anames1[1] = PTable::atom_name(e2);
    anames1[1] = "H";

    ic1.reset(2,anames1,anumbers1,xyz1);
    ic1.ic_create();
    //ic1.print_xyz();
    //ic1.print_bonds();
 
    e2 = anumbers1[1] = -1;
    anames1[1] = "X";
    ic1.reset(2,anames1,anumbers1,xyz1);
    ic1.coordn[1] = -1;
 
    int ch = climit_h[e1];
    int cl = climit_l[e1];
    int extent = ch - cl;
    for (int c1x=-1;c1x<=extent;c1x++)
    {
      int c1 = c1x + cl;
      if (c1x==-1) c1 = -1;

      printf("\n doing c1: %i for element %s \n",c1,anames1[0].c_str());
      ic1.coordn[0] = c1;

      for (int j=0;j<p[wg]+1;j++) X1[j] = X2[j] = X3[j] = 0.;
      X1[p[wg]] = X2[p[wg]] = X3[p[wg]] = 1.0;

      double norm1 = get_rowX_norm(wg,X1,ic1,1,add1,0,NULL);
      double norm2 = get_rowX_norm(wg,X2,ic1,0,NULL,1,add1);
      double norm3 = get_rowX_norm(wg,X3,ic1,1,add1,1,brk1);
      for (int k=0;k<p[wg];k++)   X3[k] -= X1[k] + X2[k];
      double norm3x = norm(X3,p[wg]+1);
      for (int k=0;k<p[wg]+1;k++) X1[k]  = X1[k] / norm1;
      for (int k=0;k<p[wg]+1;k++) X2[k]  = X2[k] / norm2;
      for (int k=0;k<p[wg]+1;k++) X3[k]  = X3[k] / norm3x;

      double val1 = knnrs[wg].predict_point(X1,KNNR_K);
      double val2 = knnrs[wg].predict_point(X2,KNNR_K);
      double val3 = knnrs[wg].predict_point(X3,KNNR_K);
      printf("  add: %4.3f brk: %4.3f add/brk: %4.3f \n",val1,val2,val3);

#if 1
      printf("  %s-%s add:",anames1[0].c_str(),anames1[1].c_str());
//      for (int k=0;k<p[wg]+1;k++)
//        printf(" %3.1f",X1[k]*10.);
      printf("    n: %4.2f \n",norm1);

      printf("  %s-%s brk:",anames1[0].c_str(),anames1[1].c_str());
//      for (int k=0;k<p[wg]+1;k++)
//        printf(" %3.1f",X2[k]*10.);
      printf("    n: %4.2f \n",norm2);

      printf("  %s(-%s) swp:",anames1[0].c_str(),anames1[1].c_str());
//      for (int k=0;k<p[wg]+1;k++)
//        printf(" %3.1f",X3[k]*10.);
      printf("  n: %4.2f \n",norm3x);
#endif
    } //loop c1 over climits
  } //loop i,j over elements

  printf("\n");

  ic1.freemem();

  delete [] qrc1;
  delete [] add1;
  delete [] brk1;
  delete [] X1;
  delete [] X2;
  delete [] X3;
#endif


  double* weight = new double[nelemf];
  double* weightp = new double[nelemf*nelemf];
  for (int i=0;i<nelemf;i++) weight[i] = 0.;
  for (int i=0;i<nelemf*nelemf;i++) weightp[i] = 0.;

  for (int i=0;i<N;i++)
  {
    int wr = i;
    int natoms1 = natoms[wr];
    string* anames2 = anames[wr];
    int* anumbers2 = anumbers[wr];
    double* xyzr2 = xyzr[wr];
    int nadd2 = nadd[wr];
    int* add2 = add[wr];
    int nbrks2 = nbrks[wr];
    int* brks2 = brks[wr];

    for (int j=0;j<nadd2;j++)
    {
      int a1 = add2[2*j+0];
      int a2 = add2[2*j+1];
      int e1 = anumbers2[a1];
      int e2 = anumbers2[a2];
      int er1 = remap[e1];
      int er2 = remap[e2];

      weight[er1] += Pr[i];
      weight[er2] += Pr[i];
      weightp[er1*nelemf+er2] += Pr[i];
      if (er1!=er2)
        weightp[er2*nelemf+er1] += Pr[i];
    }
    for (int j=0;j<nbrks2;j++)
    {
      int b1 = brks2[2*j+0];
      int b2 = brks2[2*j+1];
      int e1 = anumbers2[b1];
      int e2 = anumbers2[b2];
      int er1 = remap[e1];
      int er2 = remap[e2];

      weight[er1] += Pr[i];
      weight[er2] += Pr[i];
      weightp[er1*nelemf+er2] += Pr[i];
      if (er1!=er2)
        weightp[er2*nelemf+er1] += Pr[i];
    }
  }

  for (int i=0;i<nelemf;i++)
    weight[i] = weight[i] / N;
  for (int i=0;i<nelemf*nelemf;i++)
    weightp[i] = weightp[i] / N;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    string anam = PTable::atom_name(e1);

    printf(" element %s has weight: %4.3f \n",anam.c_str(),weight[i]);
  }
  printf("\n pair weights: \n ");
  for (int i=0;i<nelemf;i++)
    printf(" %5s",PTable::atom_name(emap[i]).c_str());
  printf("\n");
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    string anam1 = PTable::atom_name(e1);
    printf(" %2s:",anam1.c_str());
    for (int j=0;j<nelemf;j++)
    {
      //int e2 = emap[j];
      //string anam2 = PTable::atom_name(e2);

      printf(" %4.3f",weightp[i*nelemf+j]);
    }
    printf("\n");
  }


  if (allow_add_brk_elem!=NULL) delete [] allow_add_brk_elem;
  allow_add_brk_elem = new double[nelemf*nelemf];
  for (int i=0;i<nelemf*nelemf;i++) allow_add_brk_elem[i] = 0.;

  for (int i=0;i<nelemf;i++)
  for (int j=0;j<nelemf;j++)
  if (weightp[i*nelemf+j]>0.001) //CPMZ parameter
    allow_add_brk_elem[i*nelemf+j] = 1.0;

  printf(" allow_add_brk_elem: \n");
  for (int i=0;i<nelemf;i++)
  {
    for (int j=0;j<=i;j++)
      printf(" %1.0f",allow_add_brk_elem[i*nelemf+j]);
    printf("\n");
  }

  delete [] weight;
  delete [] weightp;

  return;
}



void RXNDB::ga_run(int ncycles)
{
  int stable_criterion = 3;
  double conv_crit = 0.0001;
  double lambda = LAMBDA; //was using 0.002 then 0.0005


  double mut_prob = 0.5;
  int nmut = p[0] / 6; // # of mutation attempts
//  int nmut = p[0] / 4; // # of mutation attempts
  if (nmut<2) nmut = 2;
  int nbestg = ng / 2; //keep 50%
//  int nbestg = ng / 4; //keep 25%
  if (nbestg%2==1) nbestg--;
  if (nbestg<2) nbestg = 2;
  if (ng<nbestg) nbestg = ng;


  int* bestg = new int[nbestg];
  for (int i=0;i<nbestg;i++)
    bestg[i] = i;
  double* biasl = new double[ng];

  double rmsep = 100.;
  int stable = 0;

  printf("\n\n ******** Start of GA Run! ****** \n\n");
  printf("  will perform up to %i mutations on each of %i vectors \n",nmut,ng);
  for (int n=0;n<ncycles;n++)
  {
    printf("\n  selecting the strongest, splicing and mutating \n");

    for (int i=0;i<ng;i++)
      biasl[i] = lambda * count_parameters(i);
    for (int i=0;i<ng;i++) 
      ga_error[i] += biasl[i];

#if 1
    int ndeadtot = 0;
    for (int i=0;i<ng;i++)
      ndeadtot += replace_dead_features(i);
    printf("  ndeadtot: %i \n",ndeadtot);
#endif


    selection(nbestg,bestg);

#if 0
    printf("  current best feature set: \n");
    int wgb = bestg[nbestg-1];
    for (int j=0;j<p[wgb];j++)
      rxnftrs[wgb][j].print_atoms();
#endif


    double rmse = 0.;
    int nskip = 0;
    for (int i=0;i<nbestg;i++)
    {
      int wg = bestg[i];
      if (ga_error[wg]<99999.)
        rmse += ga_error[wg]*ga_error[wg];
      else
        nskip++;
    }
    if (nbestg-nskip<1)
    {
      printf("\n  ERROR: something is wrong with regressions \n");
      printf("   reseting all features to climits \n");
      //printf(" nbestg: %i nskip: %i \n",nbestg,nskip);

      p[0] = create_features_climit(0); 
      for (int i=1;i<ng;i++)
        copy_features(i,0);
    }
    rmse = sqrt( rmse / (nbestg-nskip) );
    printf("  RMS error over best (%i): %5.4f best(#%i): %5.4f \n",nbestg-nskip,rmse,bestg[nbestg-1],ga_error[bestg[nbestg-1]]);

    printf("  GA errors (%2i):",n);
    for (int i=0;i<ng;i++)
    if (ga_error[i]<99999.)
      printf(" %4.2f",ga_error[i]);
    else
      printf("  N/A");
    printf("\n");
    printf("   penalty (x10):");
    for (int i=0;i<ng;i++)
      printf(" %4.2f",biasl[i]*10);
    printf("\n");


   //convergence criteria
    if (close_val(rmse,rmsep,conv_crit))
      stable++;
    else
      stable = 0;
    rmsep = rmse;
    if (stable>stable_criterion) break;


#if 1
    if (ncycles>1)
    {
      spawn_and_splice(nbestg,bestg);
 
      for (int i=0;i<ng;i++)
      for (int j=0;j<nmut;j++)
      {
        if (randomf(1.0)>mut_prob)
          mutate(i);
      }
    }
    else
      printf("  skipped generational update! \n");
#endif


    quiet = 1;
    //quiet = 0;
#if 1
    for (int i=0;i<ng;i++)
      ga_error[i] = create_regression(i);
#else
    for (int i=0;i<ng;i++)
      ga_error[i] = randomf(1.0);
    sleep(1.15);
#endif
    quiet = 0;

  } //loop n over ncycles


  printf("\n\n ********** GA run over ********* \n\n");

  selection(nbestg,bestg);
  for (int i=0;i<nbestg;i++)
    best_fv[i] = bestg[nbestg-1-i];
  nbest_fv = nbestg;

#if 0
  printf(" printing the best \n\n");
  for (int i=0;i<nbestg;i++)
  {
    int wg = bestg[i];
    int natoms0 = p[wg];
    printf("\n feature vector %2i with %i features \n",wg,natoms0);
    for (int j=0;j<natoms0;j++)
      rxnftrs[wg][j].print_atoms();
  }
#elif 1
  printf(" printing the best (%i) \n\n",best_fv[0]);
  int wg = best_fv[0];
  for (int j=0;j<p[wg];j++)
    rxnftrs[wg][j].print_atoms();
  ga_error[wg] = create_regression(wg);
  int nparam = 0;
  for (int j=0;j<p[wg];j++)
    nparam += rxnftrs[wg][j].count_parameters();
  printf(" feature vector has %2i total parameters \n",nparam);
#endif

  delete [] bestg;
  delete [] biasl;

  return;
}


void RXNDB::spawn_and_splice(int nbestg, int* bestg)
{
  //printf(" in spawn_and_splice for ng: %i nbestg: %i \n",ng,nbestg);

  int nneeded = ng - nbestg;
  int nmake = nneeded/2;
  int fn = nbestg;
  int on = fn + nmake;
  //printf(" ng: %i fn: %i on: %i \n",ng,fn,on);


 //first copy fittest to first positions
  int* sizes = new int[nbestg];
  RXNFTR** tmpftr = new RXNFTR*[nbestg];

  for (int i=0;i<nbestg;i++)
  {
    int wg = bestg[i];
    sizes[i] = p[wg];
    tmpftr[i] = NULL;
    copy_features(tmpftr[i],rxnftrs[wg],0,sizes[i]);
  }
  for (int i=0;i<nbestg;i++)
    copy_features(rxnftrs[i],tmpftr[i],p[i],sizes[i]);
  for (int i=0;i<nbestg;i++)
    p[i] = sizes[i];

  delete [] sizes;
  for (int i=0;i<nbestg;i++)
    delete [] tmpftr[i];
  delete [] tmpftr;

  //choose pairs at random
  for (int i=0;i<nmake;i++)
  {
    int g1 = randomi(nbestg);
    int g2 = g1;
    while (g1==g2)
      g2 = randomi(nbestg);

    int maxl = min(p[g1],p[g2]);
    int cut1 = randomi(maxl-2)+1;
    //printf(" g1,g2: %i %i cut point: %i \n",g1,g2,cut1);

    p[fn+i] = p[g2];
    p[on+i] = p[g1];

    //printf(" features new: %i %i \n",p[g2],p[g1]);

#if 0
    int natomst1 = 0;
    int natomst2 = 0;
    for (int j=0;j<p[g1];j++)
      natomst1 += rxnftrs[g1][j].natoms;
    for (int j=0;j<p[g2];j++)
      natomst2 += rxnftrs[g2][j].natoms;
    printf(" total atoms before: %i %i \n",natomst1,natomst2);
#endif

    for (int j=0;j<cut1;j++)
    {
      copy_atom(rxnftrs[fn+i][j],rxnftrs[g1][j]);
      copy_atom(rxnftrs[on+i][j],rxnftrs[g2][j]);
    }
    for (int j=cut1;j<p[g1];j++)
      copy_atom(rxnftrs[on+i][j],rxnftrs[g1][j]);
    for (int j=cut1;j<p[g2];j++)
      copy_atom(rxnftrs[fn+i][j],rxnftrs[g2][j]);

#if 0
    if (i==0)
    {
      printf("\n old features: \n");
      for (int j=0;j<p[g1];j++)
        rxnftrs[g1][j].print_atoms();
      for (int j=0;j<p[g2];j++)
        rxnftrs[g2][j].print_atoms();
      printf("\n new features: \n");
      for (int j=0;j<p[fn+i];j++)
        rxnftrs[fn+i][j].print_atoms();
      for (int j=0;j<p[on+i];j++)
        rxnftrs[on+i][j].print_atoms();
    }
#endif
#if 0
    natomst1 = 0;
    natomst2 = 0;
    for (int j=0;j<p[fn+i];j++)
      natomst1 += rxnftrs[fn+i][j].natoms;
    for (int j=0;j<p[on+i];j++)
      natomst2 += rxnftrs[on+i][j].natoms;
    printf(" total atoms after: %i %i \n",natomst1,natomst2);
#endif

  } //loop i over nmake (2 at a time)

  return;
}


void RXNDB::selection(int nbestg, int* bestg)
{
  //printf(" in selection: %i \n",nbestg);

  int* order = new int[ng];
  double* oerr = new double[ng];

  order_array(order,oerr,ga_error,ng);

  for (int i=ng-nbestg;i<ng;i++)
    bestg[i+nbestg-ng] = order[i];

#if 0
  printf("   best:");
  for (int i=0;i<nbestg;i++)
    printf(" %2i/%6.5f",bestg[i],ga_error[bestg[i]]);
  printf("\n");
#endif

  delete [] order;
  delete [] oerr;

  return;
}


void RXNDB::init_ga(int ng0)
{
  if (N<1) 
  {
    printf(" cannot init_ga_features, no data \n");
    exit(-1);
  }

  ng = ng0;
  if (ng>maxg)
  {
    printf(" ng>maxg in init_ga_features \n");
    exit(1);
  }

  create_sigmoid();
  nelemf = count_elements();

  ga_inited = 1;

  return;
}

void RXNDB::init_ga_2()
{

  printf(" initial p:");
  for (int i=0;i<ng;i++)
    printf(" %i",p[i]);
  printf("\n");

  quiet = 1;
  //quiet = 0;
  for (int i=0;i<ng;i++)
    ga_error[i] = create_regression(i);
  quiet = 0;

  printf(" GA errors:");
  for (int i=0;i<ng;i++)
    printf(" %4.2f",ga_error[i]);
  printf("\n");

  ga_inited = 2;

  return;
}
 
void RXNDB::init_ga_features(int ng0)
{
  init_ga(ng0);

#if 1
//  p[0] = create_features_climit(0);
//  p[0] = create_features_climit_ab(0);
//  p[0] = create_features_climit_abs(0);
  p[0] = create_features_climit_abs_att(0);
//  p[0] = create_features_climit_abs_att_2(0);
  printf(" ng: %i \n",ng);
  for (int i=1;i<ng;i++)
    copy_features(i,0);
#else
  for (int i=0;i<ng;i++)
    p[i] = create_features_from_db(i);
#endif

#if 0
  //random modifications
  int nmut = p[0] / 2;
  if (nmut<1) nmut = 1;
  printf("\n now performing %i mutations on each of %i vectors \n",nmut,ng);
  for (int i=0;i<ng;i++)
  for (int j=0;j<nmut;j++)
    mutate(i);
#endif

#if 0
  for (int i=0;i<ng;i++)
  {
    printf(" initial feature vector %i: \n",i);
    for (int j=0;j<p[i];j++)
      rxnftrs[i][j].print_atoms();
  }
#endif

  init_ga_2();

  return;
}

//mutates on feature vector
int RXNDB::mutate(int wg)
{
  if (rxnftrs[wg]==NULL)
  {
    printf(" rxnftrs[%i] not initialized \n",wg);
    return 0;
  }

  int nmod = 1;
  //double r1 = randomf(1.0);
  //printf(" r1: %6.4f \n",r1);
 
  int r1 = randomi(p[wg]); //feature to mod
  //printf(" r1: %i",r1);
  //if (r1>p[wg]) r1 = p[wg];

  mutate_one(wg,r1);

  return nmod;
}


//mutates one feature
void RXNDB::mutate_one(int wg, int r1)
{
  if (rxnftrs[wg]==NULL)
    return;
  if (rxnftrs[wg][r1].natoms<1)
    return;

  int natoms0 = rxnftrs[wg][r1].natoms;
  int r2 = randomi(natoms0); // atom to mod
  //printf(" r2: %i \n",r2);

//CPMZ note this
  int r3 = randomi(7); //coordn,q,abs, attached/attaching elements
//  int r3 = randomi(5); //elem,coordn,abs (swp off), atom add/delete, attached/attaching elements

#if 1
  if (r3==0)
  {
    int c1 = rxnftrs[wg][r1].get_coordn(r2);
    int r4 = randomi(2);
    int shift = 1; if (r4==0) shift = -1;
    int c2 = c1+shift;
    if (c2<-1) c2=c1-shift;
    //could prevent high coordn also
    //printf(" mutating coordn %i->%i \n",c1,c2);
    rxnftrs[wg][r1].set_coordn(r2,c2);
  }
  else if (r3==1)
  {
    int r4 = randomi(4); //25% q disable rate
    if (r4==0)
    {
      rxnftrs[wg][r1].set_qsig(r2,1000.); //disable q term
    }
    else
    {
      int r5 = randomi(2);
      int shift = 1; if (r5==0) shift = -1;
      double q1 = rxnftrs[wg][r1].get_q(r2);
      double q2 = q1 + shift * 0.2; //parameter!
      //printf(" mutating q %4.2f->%4.2f \n",q1,q2);
      rxnftrs[wg][r1].set_q(r2,q2);
      rxnftrs[wg][r1].set_qsig(r2,0.4); //parameter!
    }
  } //q mutation
  else if (1<r3 && r3<5) 
  {
    int add1,brk1,swp1;
    rxnftrs[wg][r1].get_abs(r2,add1,brk1,swp1);
    int r4 = randomi(2);
    int shift = 1; if (r4==0) shift = -1;
    //printf(" mutating abs: %2i %2i %2i -->",add1,brk1,swp1);
    if (r3==2)
    {
      add1 = add1+shift;
      if (add1<-1) add1 = add1-shift;
      if (add1> 1) add1 = add1-shift;
    }
    else if (r3==3)
    {
      brk1 = brk1+shift;
      if (brk1<-1) brk1 = brk1-shift;
      if (brk1> 1) brk1 = brk1-shift;
    }
    else if (r3==4)
    {
      swp1 = swp1+shift;
      if (swp1<-1) swp1 = swp1-shift;
      if (swp1> 1) swp1 = swp1-shift;
    }
    //printf(" %2i %2i %2i \n",add1,brk1,swp1);
    rxnftrs[wg][r1].set_abs(r2,add1,brk1,swp1);
  }
  else if (r3>=5)
  {
    int r4 = randomi(2);
    int r5 = randomi(nelemf);
    int e1 = emap[r5];
    if (r4==0)
    {
      int addAttached = 1;
      int r6 = randomi(2);
      if (rxnftrs[wg][r1].get_nattached(r2)>0 && r6==0)
          addAttached = 0;
      //else if (r6>1)
      //    addAttached = 2;
      //printf(" addAttached: %i e1: %i \n",addAttached,e1);

      if (addAttached==1)
      {
        rxnftrs[wg][r1].set_attached(r2,e1);
      }
      else if (addAttached==2)
      {
#if 1
        double q1 = randomf(1.0)-0.5;
        rxnftrs[wg][r1].set_attached(r2,e1,-1,q1);
#else
        int ch = climit_h[e1];
        int r7 = randomi(ch+1)-1;
        int c1 = r7;
        rxnftrs[wg][r1].set_attached(r2,e1,c1);
#endif
      }
      else
        rxnftrs[wg][r1].delete_attached(r2);
    } // if doing attached
    else if (r4==1)
    {
      int addAttaching = 1;
      int r6 = randomi(2);
      if (rxnftrs[wg][r2].get_nattaching(r2)>0 && r6==0)
        addAttaching = 0;
      //else if (r6>1)
      //  addAttaching = 2;
      //printf(" addAttaching: %i e1: %i \n",addAttaching,e1);

      if (addAttaching==1)
      {
        rxnftrs[wg][r1].set_attaching(r2,e1);
      }
      else if (addAttaching==2)
      {
#if 1
        double q1 = randomf(1.0)-0.5;
        printf("  mutating attaching, e1: %i q1: %4.2f \n",-1,q1);
        rxnftrs[wg][r1].set_attaching(r2,-1,-1,q1);
#else
        int ch = climit_h[e1];
        int r7 = randomi(ch+1)-1;
        int c1 = r7;
        rxnftrs[wg][r1].set_attaching(r2,e1,c1);
#endif
      }
      else
        rxnftrs[wg][r1].delete_attaching(r2);
    } // if doing attached
  }
  else
    printf("\n DEBUG: shouldn't be here \n");

#else
  if (r3==0)
  {
    int e1 = rxnftrs[wg][r1].get_element(r2);
    int e2 = e1;
    while (e2==e1)
    {
      int r4 = randomi(nelemf);
      e2 = emap[r4];
    }
    //printf(" mutating element %i->%i \n",e1,e2);
    rxnftrs[wg][r1].set_element(r2,e2);
  }
  else if (r3==1)
  {
    int c1 = rxnftrs[wg][r1].get_coordn(r2);
    int r4 = randomi(2);
    int shift = 1; if (r4==0) shift = -1;
    int c2 = c1+shift;
    if (c2<0) c2=c1-shift;
    //printf(" mutating coordn %i->%i \n",c1,c2);
    rxnftrs[wg][r1].set_coordn(r2,c2);
  }
  else if (1<r3 && r3<5) 
  {
    int add1,brk1,swp1;
    rxnftrs[wg][r1].get_abs(r2,add1,brk1,swp1);
    int r4 = randomi(2);
    int shift = 1; if (r4==0) shift = -1;
    //printf(" mutating abs: %2i %2i %2i -->",add1,brk1,swp1);
    if (r3==2)
    {
      add1 = add1+shift;
      if (add1<-1) add1 = add1-shift;
      if (add1> 1) add1 = add1-shift;
    }
    else if (r3==3)
    {
      brk1 = brk1+shift;
      if (brk1<-1) brk1 = brk1-shift;
      if (brk1> 1) brk1 = brk1-shift;
    }
    else if (r3==4)
    {
      swp1 = swp1+shift;
      if (swp1<-1) swp1 = swp1-shift;
      if (swp1> 1) swp1 = swp1-shift;
    }
    //printf(" %2i %2i %2i \n",add1,brk1,swp1);
    rxnftrs[wg][r1].set_abs(r2,add1,brk1,swp1);
  }
//CPMZ fix -- not sure what/if wrong
  else if (r3==5)
  {
    int natoms0 = rxnftrs[wg][r1].natoms;
    int r4 = randomi(2);
    int shift = 1; if (r4==0) shift = -1;
    int natoms1 = natoms0+shift;
    if (natoms1<1) natoms1 = natoms0-shift;
    if (natoms1>8) natoms1 = natoms0-shift;
    int r5 = randomi(nelemf);
    int e1 = emap[r5];
    //printf(" atom mod: %i->%i %i \n",natoms0,natoms1,e1);
    rxnftrs[wg][r1].update(natoms1);
    rxnftrs[wg][r1].set_element(-1,e1);
    //rxnftrs[wg][r1].print_atoms();
  }
  else if (r3==6)
  {
    int r4 = randomi(2);
    int r5 = randomi(nelemf);
    int e1 = emap[r5];
    if (r4==0)
    {
      int addAttached = 1;
      if (rxnftrs[wg][r1].get_nattached(r2)>0)
        addAttached = 0;
      //printf(" addAttached: %i e1: %i \n",addAttached,e1);

      if (!addAttached)
        rxnftrs[wg][r1].delete_attached(r2);
      else
      {
        //Atom at1d;
        //at1d.init();
        //at1d.set_element(e1);
        //rxnftrs[wg][r1].set_attached(r2,at1d);
        //at1d.freemem();
        rxnftrs[wg][r1].set_attached(r2,e1);
      }
    } // if doing attached
    else if (r4==1)
    {
      int addAttaching = 1;
      if (rxnftrs[wg][r1].get_nattaching(r2)>0)
        addAttaching = 0;
      //printf(" addAttaching: %i e1: %i \n",addAttaching,e1);

      if (!addAttaching)
        rxnftrs[wg][r1].delete_attaching(r2);
      else
      {
        //Atom at1d;
        //at1d.init();
        //at1d.set_element(e1);
        //rxnftrs[wg][r1].set_attached(r2,at1d);
        //at1d.freemem();
        rxnftrs[wg][r1].set_attaching(r2,e1);
      }
    } // if doing attached
  }
  else
    printf("\n DEBUG: shouldn't be here \n");
#endif


  return;
}


int RXNDB::replace_dead_features(int wg)
{
  //printf("\n feature set %i \n",wg);

#if !USE_KNNR
  double* params = lsqs[wg].params;
#else
  double* params = knnrs[wg].active;
#endif
  if (params==NULL)
    return 0;


  double THRESH = 0.00001;
  int nmut = 2;
  int ndead = 0;

#if 1
 //count ndead and return
  for (int i=0;i<p[wg];i++)
  {
    if (fabs(params[i])<THRESH)
      ndead++;
  }
  return ndead;
#endif

#if 0
  for (int i=0;i<p[wg];i++)
    rxnftrs[wg][i].print_atoms();
  printf("\n");
#endif

#if 0
  printf(" params:");
  for (int i=0;i<p[wg]+plus_b;i++)
    printf(" %5.4f",lsqs[wg].params[i]);
  printf("\n");
#endif

  for (int i=0;i<p[wg];i++)
  {
    if (fabs(params[i])<THRESH)
    {
      //printf("  feature %i is dead: %6.5f THRESH: %6.5f \n",i,lsqs[wg].params[i],THRESH);
      //printf("   feature %i is dead \n",i);
      //rxnftrs[wg][i].print_atoms();

      rxnftrs[wg][i].reset(2);

#if 1
      int ws1 = randomi(N);
      int ab1 = 0; int ab2 = 0;
      int isAdd = 0;
      if (nadd[ws1]>0 && !(randomi(2) && nbrks[ws1]))
      {
        int wa1 = randomi(nadd[ws1]);
        ab1 = add[ws1][2*wa1+0];
        ab2 = add[ws1][2*wa1+1];
        isAdd = 1;
      }
      else if (nbrks[ws1]>0)
      {
        int wb1 = randomi(nbrks[ws1]);
        ab1 = brks[ws1][2*wb1+0];
        ab2 = brks[ws1][2*wb1+1];
      }
      else
        printf(" something is off here: nadd/nbrks both zero \n");

      int e1 = anumbers[ws1][ab1];
      int e2 = anumbers[ws1][ab2];
      int c1 = coordnr[ws1][ab1];
      int c2 = coordnr[ws1][ab2];

      //printf("   ab1/2: %i %i e1/e2: %i %i c1/c2: %i %i \n",ab1,ab2,e1,e2,c1,c2);
      rxnftrs[wg][i].set_element(0,e1);
      rxnftrs[wg][i].set_element(1,e2);
      rxnftrs[wg][i].set_coordn(0,c1);
      rxnftrs[wg][i].set_coordn(1,c2);

      int add1,brk1,swp1;
      add1 = brk1 = swp1 = -1;
      if (isAdd) add1 = 1;
      else brk1 = 1;
      rxnftrs[wg][i].set_abs(0,add1,brk1,swp1);
      rxnftrs[wg][i].set_abs(1,add1,brk1,swp1);

#endif
#if 0
     //get random atoms of random structures
      int ws1 = randomi(N);
      int ws2 = randomi(N);
      int a1 = 0; int a2 = 0;
      int natoms1 = 0; int natoms2 = 0;
      while (natoms1==0)
      {
        natoms1 = natoms[ws1];
        a1 = randomi(natoms1);
      }
      int e1 = anumbers[ws1][a1];
      int e2 = e1;
      while (natoms2==0 || e2==e1)
      {
        natoms2 = natoms[ws2];
        a2 = randomi(natoms2);
        e2 = anumbers[ws2][a2];
      }

      //printf("   ws1: %i natoms1: %i a1: %i e1: %i \n",ws1,natoms1,a1,e1);
      rxnftrs[wg][i].set_element(0,e1);
      rxnftrs[wg][i].set_element(1,e2);

      int extent1 = climit_h[e1] - climit_l[e1] + 1;
      int extent2 = climit_h[e2] - climit_l[e2] + 1;
      int c1 = randomi(extent1+1)-1;
      int c2 = randomi(extent2+1)-1;
      if (c1==0) c1 = climit_h[e1];
      if (c2==0) c2 = climit_h[e2];
      rxnftrs[wg][i].set_coordn(0,c1);
      rxnftrs[wg][i].set_coordn(1,c2);
#endif
#if 1
      for (int j=0;j<nmut;j++)
        mutate_one(wg,i);
#endif

//      printf("  replaced by: \n");
//      rxnftrs[wg][i].print_atoms();
   
      ndead++;
    }
  }

  //printf("  ndead: %2i \n",ndead);

  return ndead;
}



double RXNDB::create_regression_nbo(int pr)
{
  if (Pr==NULL) create_sigmoid();

  int dim = 40; //hardcoded for now
  double* X = new double[dim*N];
  for (int i=0;i<dim*N;i++) X[i] = 0;
  //for (int i=0;i<N;i++) X[i*dim+dim-1] = 1.;

  int wg = 0;
  for (int i=0;i<N;i++)
    get_rowX_nbo(wg,i,&X[i*dim],NULL,1);

  knnrs[wg].quiet = quiet;
  knnrs[wg].ids = ids;

  printf("  using SHIFT on Pr for regression (+%3.1f kcal/mol) \n",TRAINSHIFT);
  knnrs[wg].load_values(N,dim,X,PrT2);
  knnrs[wg].load_values_print(N,Pr);
  knnrs[wg].udata = unique;
  printf("  setting unique data points \n");

  printf("  using pthresh: %5.2f \n",pthresh[wg]);
  knnrs[wg].pthresh = pthresh[wg];
  knnrs[wg].LOW_WEIGHT = 0.2;
  double error = knnrs[wg].test_points(KNNR_K);

  //printf("\n\n  note: below, the first line is actual value, subsequent are training values \n");
  int ndistant = 0;
  string scpstr = "scp";
  if (pr>1)
  {
    printf("\n printing distant \n");
    for (int i=0;i<N;i++)
    if (unique[i])
    {
      double val = knnrs[wg].errlist[i] + PrT2[i];
      if (val==1.2)
      {
        printf("  pt: %4i id: %4i val: %3.2f actual: %3.2f Ea: %5.1f \n",i,ids[i],val,Pr[i],Ea[i]);
        for (int j=0;j<dim;j++)
          printf(" %4.1f",X[i*dim+j]*10.);
        printf("\n");
        for (int l=0;l<KNNR_K;l++)
        {
          int index = knnrs[wg].knnlist[i*KNNR_K+l];
          for (int j=0;j<dim;j++)
            printf(" %4.1f",X[index*dim+j]*10.);
          printf("   %4.3f %4.3f  id: %5i/%5i",Pr[index],PrT2[index],index,ids[index]);
          printf("  dist: %4.3f \n",knnrs[wg].get_distance(&X[i*dim],&X[index*dim]));
        }
        string nstr = StringTools::int2str(ids[i],4,"0");
        scpstr += " stringfile.xyz"+nstr;
        ndistant++;
      }
    }
    printf(" %s $I \n",scpstr.c_str());
  }
  if (pr>3)
  {
    printf("\n printing positives \n");
    for (int i=0;i<N;i++)
    if (unique[i])
    {
      double val = knnrs[wg].errlist[i] + PrT2[i];
      if (val>pthresh[wg] && Pr[i]>pthresh[wg] && val<1.1)
      {
        printf("  pt: %4i id: %4i val: %3.2f actual: %3.2f Ea: %5.1f \n",i,ids[i],val,Pr[i],Ea[i]);
        for (int j=0;j<dim;j++)
          printf(" %4.1f",X[i*dim+j]*10.);
        printf("\n");
        for (int l=0;l<KNNR_K;l++)
        {
          int index = knnrs[wg].knnlist[i*KNNR_K+l];
          for (int j=0;j<dim;j++)
            printf(" %4.1f",X[index*dim+j]*10.);
          printf("   %4.3f %4.3f  id: %5i/%5i",Pr[index],PrT2[index],index,ids[index]);
          printf("  dist: %4.3f \n",knnrs[wg].get_distance(&X[i*dim],&X[index*dim]));
        }
      }
    }
  }
  if (pr>0)
  {
    scpstr = "scp";
    printf("\n printing false negatives \n");
    int nfalseneg = 0;
    for (int i=0;i<N;i++)
    if (unique[i])
    {
      double val = knnrs[wg].errlist[i] + PrT2[i];
      if (val<Pr[i]/2.0 && Pr[i]>pthresh[wg])
      {
        printf("  pt: %4i id: %4i val: %3.2f actual: %3.2f Ea: %5.1f \n",i,ids[i],val,Pr[i],Ea[i]);
        for (int j=0;j<dim;j++)
          printf(" %4.1f",X[i*dim+j]*10.);
        printf("\n");
        for (int l=0;l<KNNR_K;l++)
        {
          int index = knnrs[wg].knnlist[i*KNNR_K+l];
          for (int j=0;j<dim;j++)
            printf(" %4.1f",X[index*dim+j]*10.);
          printf("   %4.3f %4.3f  id: %5i/%5i",Pr[index],PrT2[index],index,ids[index]);
          printf("  dist: %4.3f \n",knnrs[wg].get_distance(&X[i*dim],&X[index*dim]));
        }
        string nstr = StringTools::int2str(ids[i],4,"0");
        scpstr += " stringfile.xyz"+nstr;
        nfalseneg++;
      }
    }
    printf("  found %3i false negatives \n",nfalseneg);
    printf(" %s $I \n",scpstr.c_str());
  }

  return error;
}




double RXNDB::create_regression(int wg)
{
  if (!quiet)
    printf("\n in create_regression for gene %2i \n",wg);

  if (Pr==NULL) create_sigmoid();

 //prepare features
  if (nelemf<1)
  {
    printf(" WARNING: features not ready, creating generic set \n");
    nelemf = count_elements();
    p[wg] = create_features_climit_abs(wg);
  }
  int ppb = p[wg] + plus_b;
  knnrs[wg].reset_active();

 //create X
  double* X = new double[ppb*N];
  for (int i=0;i<ppb*N;i++) X[i] = 0;
  for (int i=0;i<N;i++) X[i*ppb+ppb-1] = 1.;
#if 1
  for (int i=0;i<N;i++)
    get_rowX(wg,i,&X[i*ppb]);
#else
#endif

#if 0
  printf("\n printing X,y: \n");
  for (int i=0;i<N;i++)
  {
    for (int j=0;j<ppb;j++)
      printf(" %2.1f",X[i*ppb+j]);
    printf("  %4.3f   id: %3i/%3i \n",Pr[i],i,ids[i]);
  }
  printf("\n");
#endif

  //printf(" quiet? %i \n",quiet);
#if !USE_KNNR
  lsqs[wg].quiet = quiet;
  double error = lsqs[wg].do_least_squares(X,ppb,Pr,N);
#else
  knnrs[wg].quiet = quiet;
  knnrs[wg].ids = ids;

#if TEST_GA && MAX_N
  int N1 = N;
  if (N1>MAX_N) N1 = MAX_N;
  knnrs[wg].load_values(N,ppb,X,Pr);
  knnrs[wg].knn_N = N1;
#else
  knnrs[wg].load_values(N,ppb,X,Pr);
#endif

  knnrs[wg].pthresh = pthresh[wg];
  double error = knnrs[wg].test_points(KNNR_K);
#endif

#if TEST_GA && USE_KNNR && 1
  if (!quiet)
  {
    printf("\n after regression, printing positive results \n");
    for (int i=0;i<knnrs[wg].npts;i++)
    {
      double val = knnrs[wg].errlist[i] + Pr[i];
      if (val > pthresh[wg])
      {
        printf("  pt: %4i id: %4i val: %3.2f actual: %3.2f \n",i,ids[i],val,Pr[i]);
        //for (int j=0;j<nadd[i];j++)
        //  printf("   adding: %i %i \n",add[i][2*j+0]+1,add[i][2*j+1]+1);
        //for (int j=0;j<nbrks[i];j++)   
        //  printf("   breaking: %i %i \n",brks[i][2*j+0]+1,brks[i][2*j+1]+1);

        for (int j=0;j<ppb;j++)
          printf(" %2.1f",X[i*ppb+j]*10.);
        printf("  %4.3f   id: %3i/%3i \n",Pr[i],i,ids[i]);
        for (int l=0;l<KNNR_K;l++)
        {
          int index = knnrs[wg].knnlist[i*KNNR_K+l];
          for (int j=0;j<ppb;j++)
            printf(" %2.1f",X[index*ppb+j]*10.);
          printf("  %4.3f   id: %4i/%4i",Pr[index],index,ids[index]);
          printf("  dist: %4.3f \n",knnrs[wg].get_distance(&X[i*ppb],&X[index*ppb]));
        }
      } //if positive

    } //loop i over npts
    printf("\n after regression, printing false negatives \n");
    int nfalseneg = 0;
    for (int i=0;i<knnrs[wg].npts;i++)
    {
      double val = knnrs[wg].errlist[i] + Pr[i];
      if (Pr[i]>0.4 && val < Pr[i]/2.) //new false neg criterion, was 2.5
      {
        printf("  pt: %4i id: %4i val: %3.2f actual: %3.2f \n",i,ids[i],val,Pr[i]);
        for (int j=0;j<ppb;j++)
          printf(" %2.1f",X[i*ppb+j]*10.);
        printf("  %4.3f   id: %4i/%4i \n",Pr[i],i,ids[i]);
        for (int l=0;l<KNNR_K;l++)
        {
          int index = knnrs[wg].knnlist[i*KNNR_K+l];
          for (int j=0;j<ppb;j++)
            printf(" %2.1f",X[index*ppb+j]*10.);
          printf("  %4.3f   id: %3i/%3i",Pr[index],index,ids[index]);
          printf("  dist: %4.3f \n",knnrs[wg].get_distance(&X[i*ppb],&X[index*ppb]));
        }
        nfalseneg++;
      } //if false negative

    } //loop i over npts
    printf("  found %i false negatives \n",nfalseneg);
  } //if  !quiet

#endif

  if (error<99999. && !quiet)
  {
    printf("   RMS error: %4.3f \n",error);
//    printf("   total unsigned error: %4.3f \n",error);
//    printf(" average unsigned error: %4.3f \n",error/N);
  }
  else if (!quiet)
    printf("  regression failed \n");


  delete [] X;

  return error;
}


void RXNDB::create_sigmoid()
{
  double kT = temperature*(0.6/300.); //need precision 
  double Eref = 60*kT;
  double alpha_s = 10.*kT; //was 5kT

  printf(" T: %4.2f K  kT: %3.2f kcal/mol  Eref: %3.2f kcal/mol  sigmoid_alpha: %3.2f  N: %2i \n",temperature,kT,Eref,alpha_s,N);

  if (Pr!=NULL) delete [] Pr;
  Pr = new double[N];
  PrT2 = new double[N];
  for (int i=0;i<N;i++)
    Pr[i] = sigmoid(Ea[i],Eref,alpha_s);
  for (int i=0;i<N;i++)
    PrT2[i] = sigmoid(Ea[i],Eref+TRAINSHIFT,alpha_s);

  printf(" sigmoid Pr:");
  for (int i=0;i<N;i++)
    printf(" %5.3f",Pr[i]);
  printf("\n");

  int npass = 0;
  for (int i=0;i<N;i++)
  if (Pr[i]>0.5)
    npass++;
  double passrate = 100.0*npass/N;
  printf(" exact pass rate (%4i of %4i): %4.1f%% \n",npass,N,passrate);

  return;
}



void RXNDB::get_rowX_nbo(int wg, int wr, double* X, int* M, int sort)
{
 //need to finish this feature set!
 //expand to 5add 5brk

#define SPRSCALE 1.

  int dim = 40;
  for (int i=0;i<dim;i++) X[i] = 0.;

  double* abpol = new double[5];
  double* abspr = new double[10];
  double* abq = new double[10];
  int* ms = new int[10];
  for (int k=0;k<10;k++) ms[k] = 0;
  for (int k=0;k<10;k++) abq[k] = abspr[k] = 0.;

  int nbrk1 = min(nbosr[wr].nb,5);
  int nadd1 = min(nbosp[wr].nb,5);
 // if (nbosr[wr].nb>5) printf("   wr: %4i ids: %4i nbrk1: %2i Ea: %5.1f \n",wr,ids[wr],nbosr[wr].nb,Ea[wr]);
 // if (nbosp[wr].nb>5) printf("   wr: %4i ids: %4i nadd1: %2i Ea: %5.1f \n",wr,ids[wr],nbosp[wr].nb,Ea[wr]);

  //addq[8](addpol) - addspr[8] - brkq[8](brkpol) - brkspr[8]

  for (int j=0;j<nadd1;j++)
  {
    int a1  = nbosp[wr].blist[j];
    int at1 = nbosp[wr].bmo_atoms[2*a1+0]-1;
    int at2 = nbosp[wr].bmo_atoms[2*a1+1]-1;
    int tz = 0;
    if (at2<0) { at2 = at1; tz = 1; }
    //printf(" a1: %2i at1/2: %2i %2i (%2s %2s) \n",a1,at1,at2,anames[wr][at1].c_str(),anames[wr][at2].c_str());

    abq[2*j+0] = nbosr[wr].q[at1];
    abq[2*j+1] = nbosr[wr].q[at2];
    abspr[2*j+0] = nbosr[wr].spratio[at1];
    abspr[2*j+1] = nbosr[wr].spratio[at2];
    ms[2*j+0] = at1 + 1;
    ms[2*j+1] = at2 + 1;
    if (tz) abq[2*j+1] = 0.;
    if (tz) abspr[2*j+1] = 0.;
    if (tz) ms[2*j+1] = 0;
    if (sort) sort_2p(&abspr[2*j],&abq[2*j]);
//    if (sort) sort_2p(&abq[2*j],&abspr[2*j]);
  } 
  if (sort)
  {
    //printf("  before sort: %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f \n",abspr[0],abspr[1],abspr[2],abspr[3],abspr[4],abspr[5],abspr[6],abspr[7],abspr[8],abspr[9]);
    if (nadd1==2) sort_n2p(abspr,abq);
    if (nadd1==3) sort_n3p(abspr,abq);
    if (nadd1==4) sort_n4p(abspr,abq);
    if (nadd1==5) sort_n5p(abspr,abq);
    //printf("   after sort: %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f \n",abspr[0],abspr[1],abspr[2],abspr[3],abspr[4],abspr[5],abspr[6],abspr[7],abspr[8],abspr[9]);
  }

  for (int j=0;j<5;j++)
  {
    X[0+2*j+0] = abq[2*j+0];
    X[0+2*j+1] = abq[2*j+1];
    X[10+2*j+0] = abspr[2*j+0]/SPRSCALE;
    X[10+2*j+1] = abspr[2*j+1]/SPRSCALE;
    if (M!=NULL)
    {
      M[0+2*j+0] = ms[2*j+0];
      M[0+2*j+1] = ms[2*j+1];
      M[10+2*j+0] = ms[2*j+0];
      M[10+2*j+1] = ms[2*j+1];
    }
  } //first 20 values of X

  for (int j=0;j<nbrk1;j++)
  {
    int b1  = nbosr[wr].blist[j];
    int bt1 = nbosr[wr].bmo_atoms[2*b1+0]-1;
    int bt2 = nbosr[wr].bmo_atoms[2*b1+1]-1;
    int tz = 0;
    if (bt2<0) { bt2 = bt1; tz = 1; }
    //printf(" b1: %2i bt1/2: %2i %2i (%2s %2s) \n",b1,bt1,bt2,anames[wr][bt1].c_str(),anames[wr][bt2].c_str());

    abq[2*j+0] = nbosr[wr].q[bt1];
    abq[2*j+1] = nbosr[wr].q[bt2];
    abspr[2*j+0] = nbosr[wr].spratio[bt1];
    abspr[2*j+1] = nbosr[wr].spratio[bt2];
    ms[2*j+0] = bt1 + 1;
    ms[2*j+1] = bt2 + 1;
    if (tz) abq[2*j+1] = 0.;
    if (tz) abspr[2*j+1] = 0.;
    if (tz) ms[2*j+1] = 0;
    if (sort) sort_2p(&abspr[2*j],&abq[2*j]);
//    if (sort) sort_2p(&abq[2*j],&abspr[2*j]);
  }
  if (sort)
  {
    if (nbrk1==2) sort_n2p(abspr,abq);
    if (nbrk1==3) sort_n3p(abspr,abq);
    if (nbrk1==4) sort_n4p(abspr,abq);
    if (nbrk1==5) sort_n5p(abspr,abq);
  }

  for (int j=0;j<5;j++)
  {
    X[20+2*j+0] = abq[2*j+0];
    X[20+2*j+1] = abq[2*j+1];
    X[30+2*j+0] = abspr[2*j+0]/SPRSCALE;
    X[30+2*j+1] = abspr[2*j+1]/SPRSCALE;
    if (M!=NULL)
    {
      M[20+2*j+0] = ms[2*j+0];
      M[20+2*j+1] = ms[2*j+1];
      M[30+2*j+0] = ms[2*j+0];
      M[30+2*j+1] = ms[2*j+1];
    }
  }


  delete [] abpol;
  delete [] abq;
  delete [] abspr;
  delete [] ms;

#if 0
  if (sort>0)
  {
    double norm1 = norm(X,dim); //norm over +b is important
    if (norm1>0.)
    for (int i=0;i<dim;i++)
      X[i] = X[i] / norm1;
    else
    for (int i=0;i<dim;i++)
      X[i] = 0.;
  }
#endif

  return;
}

void RXNDB::get_rowX(int wg, int wr, double* X)
{
  int natoms1 = natoms[wr];
  string* anames1 = anames[wr];
  int* anumbers1 = anumbers[wr];
  double* xyzr1 = xyzr[wr];
  int nadd1 = nadd[wr];
  int* add1 = add[wr];
  int nbrks1 = nbrks[wr];
  int* brks1 = brks[wr];

  ICoord ic1;
  ic1.init(natoms1,anames1,anumbers1,xyzr1);
  ic1.q = qrc[wr];

  //printf("\n rxn %3i: \n",wr);
  get_rowX(wg,X,ic1,nadd1,add1,nbrks1,brks1);

  ic1.freemem();

  return;
}

void RXNDB::get_rowX(int wg, double* X, ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1)
{
  int size = p[wg];
  int sizepb = size + plus_b;

  int check_active = 0;
  if (knnrs[wg].active!=NULL) check_active = 1;

 //note: get_value_1 also available
  if (check_active)
  {
    for (int i=0;i<size;i++)
    if (knnrs[wg].active[i])
      X[i] = rxnftrs[wg][i].get_value_0(ic1,nadd1,add1,nbrks1,brks1);
  }
  else
  {
    for (int i=0;i<size;i++)
      X[i] = rxnftrs[wg][i].get_value_0(ic1,nadd1,add1,nbrks1,brks1);
  }

  //printf(" here get_rowX \n"); fflush(stdout);

#if 1
  double norm1 = norm(X,sizepb); //norm over +b is important
  if (norm1>0.)
  for (int i=0;i<sizepb;i++)
    X[i] = X[i] / norm1;
  else
  for (int i=0;i<sizepb;i++)
    X[i] = 0.;
#endif

  return;
}

double RXNDB::get_rowX_norm(int wg, double* X, ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1)
{
  int size = p[wg];
  int sizepb = size + plus_b;
#if 1
  for (int i=0;i<size;i++)
  if (knnrs[wg].active[i])
    X[i] = rxnftrs[wg][i].get_value_0(ic1,nadd1,add1,nbrks1,brks1);
#else
  for (int i=0;i<p[wg];i++)
  if (knnrs[wg].active[i])
    X[i] = rxnftrs[wg][i].get_value_1(ic1,nadd1,add1,nbrks1,brks1);
#endif

  //printf(" here get_rowX \n"); fflush(stdout);

  double norm1 = norm(X,sizepb); //norm over +b is important

  return norm1;
}


double RXNDB::eval_prob(int natoms0, string* anames0, int* anumbers0, double* xyzr0, int nadd0, int* add0, int nbrks0, int* brks0, double* qrc0)
{
  double val = 0.;

  //printf("  eval_prob nadd: %i nbrks: %i \n",nadd0,nbrks0);
  //printf("  wg: %i \n",best_fv[0]);
  //for (int i=0;i<nadd0;i++)
  //  printf("   adding: %i %i \n",add0[2*i+0]+1,add0[2*i+1]+1);
  //for (int i=0;i<nbrks0;i++)   
  //  printf("   breaking: %i %i \n",brks0[2*i+0]+1,brks0[2*i+1]+1);
  //printf(" qrc0:");
  //for (int i=0;i<natoms0;i++)
  //  printf(" %4.2f",qrc0[i]);
  //printf("\n");

  if (p[0]<1)
  {
    printf("  cannot make estimate, p=0 \n");
    return -999.;
  }

  int ppbmax = 0;
  for (int i=0;i<nbest_fv;i++)
  if (ppbmax<p[i]+1)
    ppbmax = p[i]+1;
  double* X = new double[ppbmax];
  ICoord ic1;
  ic1.init(natoms0,anames0,anumbers0,xyzr0);
  ic1.q = qrc0;

  double valt = 0;
  int nreg = 1;
  for (int i=0;i<nreg;i++)
  {
    int wg = best_fv[i];
    //printf("  wg: %2i",wg);

    int ppb = p[wg] + plus_b;
    for (int i=0;i<p[wg];i++) X[i] = 0.;
    X[ppb-1] = 1.;

    get_rowX(wg,X,ic1,nadd0,add0,nbrks0,brks0);

    double mag = norm(X,ppb);
    if (mag>0.) //make sure vector is not zero
    {
#if !USE_KNNR
      val = lsqs[wg].estimate_point(X);
#else
      val = knnrs[wg].predict_point(X,KNNR_K);
      if (knnrs[wg].sumw < 0.01)
        val = 1.25;
#endif
    }
    else
      val = 1.3;

#if USE_KNNR && 1
//CPMZ here
    printf(" val: %4.3f  (wg: %i) \n",val,wg);
    for (int j=0;j<ppb;j++)
      printf(" %2.1f",X[j]*10.);
    printf("\n");
    for (int l=0;l<KNNR_K;l++)
    {
      int index = knnrs[wg].knnlist[l];
      double* X1 = &knnrs[wg].X[index*ppb];

      for (int j=0;j<ppb;j++)
        printf(" %2.1f",X1[j]*10.);
      printf("  %4.3f   id: %3i/%3i",Pr[index],index,ids[index]);
      printf("  dist: %4.3f \n",knnrs[wg].get_distance(X,X1));
    }
#endif

    if (val>pthresh[wg])
    {
      //printf(" %4.3f > %4.3f valt++ \n",val,pthresh[wg]);
      valt++;
    }
  }

  delete [] X;
  ic1.freemem();

  valt = valt / nreg;

  return valt;
}



int RXNDB::check_duplicate(int ws0, int natoms0, int* anumbers0, double* energies,
           int nadd0, int* add0, int nbrks0, int* brks0)
{
  int found = 0;
//  return found;

  int wf = -1;
  for (int i=0;i<N;i++)
  if (!found)
  if (natoms0==natoms[i])
  if (nadd0==nadd[i])
  if (nbrks0==nbrks[i])
  {
    //same natoms, nadd, nbrks
    int anum = 0;
    for (int j=0;j<natoms0;j++)
    if (anumbers0[j]==anumbers[i][j])
      anum++;
    if (anum==natoms0)
    {
      //check that adds are the same
      int nsame = 0;
      int* same = new int[nadd0+nbrks0];
      for (int j=0;j<nadd0;j++) same[j] = 0;
      for (int j=0;j<nadd0;j++)
      {
        int a1n = add0[2*j+0];
        int a2n = add0[2*j+1];
        for (int k=0;k<nadd0;k++)
        if (!same[k])
        {
          int a1 = add[i][2*k+0];
          int a2 = add[i][2*k+1];
          if ((a1==a1n && a2==a2n) || (a1==a2n && a2==a1n))
          {
            nsame++;
            same[k] = 1;
            break;
          }
        }
      } //loop j over nadd0
      //then check brks are the same
      for (int j=0;j<nbrks0;j++) same[j] = 0;
      for (int j=0;j<nbrks0;j++)
      {
        int b1n = brks0[2*j+0];
        int b2n = brks0[2*j+1];
        for (int k=0;k<nbrks0;k++)
        if (!same[k])
        {
          int b1 = brks[i][2*k+0];
          int b2 = brks[i][2*k+1];
          if ((b1==b1n && b2==b2n) || (b1==b2n && b2==b1n))
          {
            nsame++;
            same[k] = 1;
            break;
          }
        }
      }
      delete [] same;
      if (nsame==nadd0+nbrks0)
      {
        wf = i;
        found = 1;
        break;
      }
    } //if anumbers match
  } //loop i over N points, if prelim match

  if (found)
  {
    //printf("   new ts %3i same as %3i \n",ws0,ids[wf]);
#if 0
    for (int i=0;i<nadd0;i++)
      printf(" add: %i-%i %i-%i \n",add0[2*i+0],add0[2*i+1],add[wf][2*i+0],add[wf][2*i+1]);
    for (int i=0;i<nbrks0;i++)
      printf(" brk: %i-%i %i-%i \n",brks0[2*i+0],brks0[2*i+1],brks[wf][2*i+0],brks[wf][2*i+1]);
#endif
    found = wf;
  }
  else
    found = -1;

  return found;
}


int RXNDB::add_ts(int ws, int natoms0, string* anames0, int* anumbers0, double* xyzr0, double* xyzp0, double* xyzts0, double* energies,
           int nadd0, int* add0, int nbrks0, int* brks0, double* qrc0)
{
  //printf(" in add_ts (natoms: %2i) \n",natoms0);
  //printf(" nadd: %i nbrks: %i \n",nadd0,nbrks0);
  //printf(" energies: %4.1f %4.1f %4.1f \n",energies[0],energies[2],energies[1]); //r/ts/p
  //fflush(stdout);

  lastAdd = -1;
  if (natoms0<1)
  {
    printf("  too few atoms (%i), cannot add_ts \n",natoms0);
    return 1;
  }
  if (nadd0==0 && nbrks0==0)
  {
    //printf(" no adds or breaks, not adding to db \n");
    return 1;
  }

  int isnew = 1;
  int wsn = N;

#if DO_NOT_ADD_DUPLICATES
  int dupl = check_duplicate(ws,natoms0,anumbers0,energies,nadd0,add0,nbrks0,brks0);
#else
  int dupl = -1;
#endif
  if (dupl>-1)
  {
    //printf(" Ea/Ea: %3.1f %3.1f \n",energies[2]-energies[0],Ea[dupl]);
    if (Ea[dupl]>energies[2]-energies[0])
    {
      wsn = dupl;
      isnew = 0;
    }
    else 
      return 0; //return no error
  }


  ICoord ic1;
  ic1.init(natoms0,anames0,anumbers0,xyzr0);

 //save descriptors
  if (isnew)
    alloc_db(natoms0);

  lastAdd = wsn;
  ids[wsn] = ws;
  natoms[wsn] = natoms0;
  for (int i=0;i<natoms0;i++)
    anames[wsn][i] = anames0[i]; 
  for (int i=0;i<natoms0;i++)
    anumbers[wsn][i] = anumbers0[i]; 
  for (int i=0;i<3*natoms0;i++)
    xyzr[wsn][i] = xyzr0[i];
  for (int i=0;i<3*natoms0;i++)
    xyzp[wsn][i] = xyzp0[i];
  for (int i=0;i<natoms0;i++)
    coordnr[wsn][i] = ic1.coordn[i];
  nadd[wsn] = nadd0;
  for (int i=0;i<2*nadd0;i++)
    add[wsn][i] = add0[i];
  nbrks[wsn] = nbrks0;
  for (int i=0;i<2*nbrks0;i++)
    brks[wsn][i] = brks0[i];
  Ea[wsn] = energies[2] - energies[0]; //assumes kcal/mol
  Erxn[wsn] = energies[1] - energies[0]; 

  for (int i=0;i<natoms0;i++)
    qrc[wsn][i] = qrc0[i];

#if 0
  printf("  added qrc0:");
  for (int i=0;i<natoms0;i++)
    printf(" %4.2f",qrc0[i]);
  printf("\n");
#endif

  if (isnew)
    printf("     added TS(#%i/id:%i) with Ea: %4.1f Erxn: %4.1f \n",N,ws,Ea[wsn],Erxn[wsn]);
  else
    printf("     replaced TS(#%i/id:%i) with Ea: %4.1f Erxn: %4.1f \n",wsn,ws,Ea[wsn],Erxn[wsn]);
 
  if (isnew)
    N++;

  ic1.freemem();

  return 0;
}


int RXNDB::add_ts_xyz(int ws, int natoms0, string* anames0, int* anumbers0, double* xyzr0, double* xyzp0, double* xyzts0, double* energies, double* qrc0)
{
  //printf(" in add_ts_xyz (natoms: %2i) \n",natoms0);

  if (natoms0<1)
  {
    printf("  too few atoms (%i), cannot add_ts \n",natoms0);
    return 1;
  }
  if (anames0==NULL || anumbers0==NULL || xyzr0==NULL || xyzp0==NULL || xyzts0==NULL || energies==NULL)
  {
    printf("  pointer problem in add_ts_xyz \n");
    return 1;
  }

  int error = 0;

  ICoord ic1,ic2;
  ic1.init(natoms0,anames0,anumbers0,xyzr0);
  ic2.init(natoms0,anames0,anumbers0,xyzp0);

  int nadd1 = 0;
  int* add1 = new int[20];
  int nbrks1 = 0;
  int* brks1 = new int[20];

 //search for differences in bonding
  for (int i=0;i<ic2.nbonds;i++)
  {
    int a1 = ic2.bonds[i][0];
    int a2 = ic2.bonds[i][1];
    if (!ic1.bond_exists(a1,a2))
    {
      if (nadd1>9) 
      {
        printf(" cannot add, too many bond differences \n");
        error = 2;
        break;
      }
      add1[2*nadd1+0] = a1;
      add1[2*nadd1+1] = a2;
      nadd1++;
      //printf(" new bond: %i %i \n",a1+1,a2+1);
    }
  }

  if (!error)
  for (int i=0;i<ic1.nbonds;i++)
  {
    int b1 = ic1.bonds[i][0];
    int b2 = ic1.bonds[i][1];
    if (!ic2.bond_exists(b1,b2))
    {
      if (nbrks1>9) 
      {
        printf(" cannot add, too many bond differences \n");
        error = 2;
        break;
      }
      brks1[2*nbrks1+0] = b1;
      brks1[2*nbrks1+1] = b2;
      nbrks1++;
      //printf(" broken bond: %i %i \n",b1+1,b2+1);
    }
  }

  //printf("  about to add_ts \n"); fflush(stdout);
  if (!error)
    error = add_ts(ws,natoms0,anames0,anumbers0,xyzr0,xyzp0,xyzts0,energies,nadd1,add1,nbrks1,brks1,qrc0);

  delete [] add1;
  delete [] brks1;

  ic1.freemem();
  ic2.freemem();

  return error;
}



int RXNDB::add_ts_xyz_test(int ws, int natoms0, string* anames0, int* anumbers0, double* xyzr0, double* xyzp0, double* xyzts0, double* energies, double* qrc0, int nadd0, int* add0, int nbrks0, int* brks0)
{
  //printf(" in add_ts_xyz_test (natoms: %2i) \n",natoms0);

  if (natoms0<1)
  {
    printf("  too few atoms (%i), cannot add_ts \n",natoms0);
    return 1;
  }
  if (anames0==NULL || anumbers0==NULL || xyzr0==NULL || xyzp0==NULL || xyzts0==NULL || energies==NULL)
  {
    printf("  pointer problem in add_ts_xyz \n");
    return 1;
  }

  int error = 0;

  ICoord ic1,ic2;
  ic1.init(natoms0,anames0,anumbers0,xyzr0);
  ic2.init(natoms0,anames0,anumbers0,xyzp0);

  int nadd1 = 0;
  int* add1 = new int[20];
  int nbrks1 = 0;
  int* brks1 = new int[20];

 //search for differences in bonding
  for (int i=0;i<ic2.nbonds;i++)
  {
    int a1 = ic2.bonds[i][0];
    int a2 = ic2.bonds[i][1];
    if (!ic1.bond_exists(a1,a2))
    {
      if (nadd1>9) 
      {
        printf(" cannot add, too many bond differences \n");
        error = 2;
        break;
      }
      add1[2*nadd1+0] = a1;
      add1[2*nadd1+1] = a2;
      nadd1++;
      //printf("   new bond: %i %i \n",a1+1,a2+1);
    }
  }

  if (!error)
  for (int i=0;i<ic1.nbonds;i++)
  {
    int b1 = ic1.bonds[i][0];
    int b2 = ic1.bonds[i][1];
    if (!ic2.bond_exists(b1,b2))
    {
      if (nbrks1>9) 
      {
        printf(" cannot add, too many bond differences \n");
        error = 2;
        break;
      }
      brks1[2*nbrks1+0] = b1;
      brks1[2*nbrks1+1] = b2;
      nbrks1++;
      //printf("   broken bond: %i %i \n",b1+1,b2+1);
    }
  }

  //printf("  about to add_ts \n"); fflush(stdout);
  if (!error)
    error = add_ts(ws,natoms0,anames0,anumbers0,xyzr0,xyzp0,xyzts0,energies,nadd1,add1,nbrks1,brks1,qrc0);

 //compare ISOMERS to actual changes
  int nsim = 0;
  if (nadd0==nadd1)
    nsim += nadd0;
  if (nbrks0==nbrks1)
    nsim += nbrks0;
  //printf("   same number of add/brk? %i addbrk(max): %i \n",nsim,max(nadd0+nbrks0,nadd1+nbrks1));

//CPMZ here

 //did the moves requested happen?
  int naf = 0;
  for (int i=0;i<nadd0;i++)
  for (int j=0;j<nadd1;j++)
  {
    if (add1[2*j+0]==add0[2*i+0] && add1[2*j+1]==add0[2*i+1])
      naf++;
    if (add1[2*j+0]==add0[2*i+1] && add1[2*j+1]==add0[2*i+0])
      naf++;
  }
  int nbf = 0;
  for (int i=0;i<nbrks0;i++)
  for (int j=0;j<nbrks1;j++)
  {
    if (brks1[2*j+0]==brks0[2*i+0] && brks1[2*j+1]==brks0[2*i+1])
      nbf++;
    if (brks1[2*j+0]==brks0[2*i+1] && brks1[2*j+1]==brks0[2*i+0])
      nbf++;
  }
  int sa = 0;
  int sb = 0;
  if (naf==nadd0) sa = 1;
  if (nbf==nbrks0) sb = 1;
  //printf("   same add/brk? %i %i fr: %4.3f \n",sa,sb,0.5*(sa+sb));

  if (lastAdd>-1)
    similar[lastAdd] = 0.5*(sa + sb);



  delete [] add1;
  delete [] brks1;

  ic1.freemem();
  ic2.freemem();

  return error;
}


int RXNDB::count_elements()
{
  printf("\n counting elements and creating map \n");
  int nfound = 0;

  int* exists = new int[nelem];
  for (int i=0;i<nelem;i++) exists[i] = 0;

  for (int i=0;i<N;i++)
  for (int j=0;j<natoms[i];j++)
  {
    int e1 = anumbers[i][j] - 1;
    exists[e1] = 1;
  }


 //emap is sequentially ordered elements
  if (emap==NULL) emap = new int[nelem];
  for (int i=0;i<nelem;i++) emap[i] = 0;
  for (int i=0;i<nelem;i++)
  if (exists[i])
    emap[nfound++] = i+1;

 //remap[atomicnum] = index of emap
  if (remap==NULL) remap = new int[nelem];
  for (int i=0;i<nelem;i++) remap[i] = -1;
  for (int i=0;i<nfound;i++)
    remap[emap[i]] = i;


  delete [] exists;

#if 1
  printf(" found elements:");
  for (int i=0;i<nfound;i++)
    printf(" %2i",emap[i]);
  printf("\n");
#endif
#if 0
  printf(" remap indices:");
  for (int i=0;i<nelem;i++)
    printf(" %2i",remap[i]);
  printf("\n");
#endif

  if (nfound==0) 
  { 
    printf(" ERROR: couldn't find any elements \n");
    exit(-1);
  }
  else if (nfound==1)
  {
    printf(" ERROR: may not work with only 1 element \n");
    exit(-1);
  }
  return nfound;
}


int RXNDB::count_parameters(int wg)
{
  int nparam = 0;

  for (int i=0;i<p[wg];i++)
    nparam += rxnftrs[wg][i].count_parameters();

  return nparam;
}


int RXNDB::create_features_from_db(int wg)
{
  printf("\n creating element/coordn features from db for vector %i \n",wg); 

  if (rxnftrs!=NULL)
  if (rxnftrs[wg]!=NULL) 
  {
    //printf(" freeing mem in rxnftrs \n"); fflush(stdout);
    for (int i=0;i<p[wg];i++)
      rxnftrs[wg][i].freemem();
    delete [] rxnftrs[wg];
    //printf(" done freeing mem in rxnftrs \n"); fflush(stdout);
  }

  int ntotal = 0;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    ntotal += climit_h[e1] - climit_l[e1] + 1;
  }
  //printf(" will create %i features \n",ntotal);
  rxnftrs[wg] = new RXNFTR[ntotal];

  int* structs = new int[ntotal];
  int nfound = 0;
  int nnew = 0;
  for (int i=0;i<ntotal;i++)
  {
   //locate unique structure from db
    int ws1 = -1;
    int fnew = 0;
    while (!fnew)
    {
      ws1 = randomi(N);
      fnew = 1;
      for (int j=0;j<nfound;j++)
      if (ws1==structs[j])
        fnew = 0;
    }
    structs[nfound++] = ws1;

    int size = 2*nadd[ws1] + 2*nbrks[ws1];
    printf(" size is: %i \n",size);
    if (size<1)
      printf(" something is off here: nadd/nbrks both zero \n");

  //change this to collect all relevant atoms?
    rxnftrs[wg][nnew].init(size);

    int ab1 = 0; int ab2 = 0;
    int add1,brk1,swp1;
   //loop over nadd, nbrk, create features
    for (int j=0;j<nadd[ws1];j++)
    {
      int index1 = 2*j+0;
      int index2 = 2*j+1;

      ab1 = add[ws1][index1];
      ab2 = add[ws1][index2];

      int e1 = anumbers[ws1][ab1];
      int e2 = anumbers[ws1][ab2];
      int c1 = coordnr[ws1][ab1];
      int c2 = coordnr[ws1][ab2];

      //printf("   ab1/2: %i %i e1/e2: %i %i c1/c2: %i %i \n",ab1,ab2,e1,e2,c1,c2);

      rxnftrs[wg][nnew].set_element(index1,e1);
      rxnftrs[wg][nnew].set_element(index2,e2);
      rxnftrs[wg][nnew].set_coordn(index1,c1);
      rxnftrs[wg][nnew].set_coordn(index2,c2);
 
      add1 = 1; brk1 = -1; swp1 = -1;
      rxnftrs[wg][nnew].set_abs(index1,add1,brk1,swp1);
      rxnftrs[wg][nnew].set_abs(index2,add1,brk1,swp1);
    } //add features
    for (int j=0;j<nbrks[ws1];j++)
    {
      ab1 = brks[ws1][2*j+0];
      ab2 = brks[ws1][2*j+1];

      int e1 = anumbers[ws1][ab1];
      int e2 = anumbers[ws1][ab2];
      int c1 = coordnr[ws1][ab1];
      int c2 = coordnr[ws1][ab2];
      //printf("   ab1/2: %i %i e1/e2: %i %i c1/c2: %i %i \n",ab1,ab2,e1,e2,c1,c2);

      int index1 = 2*nadd[ws1]+2*j+0;
      int index2 = 2*nadd[ws1]+2*j+1;

      rxnftrs[wg][nnew].set_element(index1,e1);
      rxnftrs[wg][nnew].set_element(index2,e2);
      rxnftrs[wg][nnew].set_coordn(index1,c1);
      rxnftrs[wg][nnew].set_coordn(index2,c2);
 
      add1 = -1; brk1 = 1; swp1 = -1;
      rxnftrs[wg][nnew].set_abs(index1,add1,brk1,swp1);
      rxnftrs[wg][nnew].set_abs(index2,add1,brk1,swp1);
    } //brk features

    nnew++;
  }

  delete [] structs;

  nfound = nnew;
  plus_b = 1;

  printf("\n found %2i features (+b: %i)\n",nfound,plus_b); 

  return nfound;
}


int RXNDB::create_features_climit(int wg)
{
  printf("\n creating element/coordn features via climits \n"); 

  int nfound = 0;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;
    nfound += extent;

    printf("  element %2i extent: %i \n",e1,extent);
  }
#if NULL_ATOMS
  int nt = NULL_ATOMS;
  nfound += nelemf*nt;
#endif

  if (rxnftrs!=NULL)
  if (rxnftrs[wg]!=NULL) 
  {
    //printf(" freeing mem in rxnftrs \n"); fflush(stdout);
    for (int i=0;i<p[wg];i++)
      rxnftrs[wg][i].freemem();
    delete [] rxnftrs[wg];
    //printf(" done freeing mem in rxnftrs \n"); fflush(stdout);
  }
  rxnftrs[wg] = new RXNFTR[nfound];


 //create feature vectors for each
  int nnew = 0;
#if NULL_ATOMS
  for (int i=0;i<nt*nelemf;i++)
  {
    int e1 = -1;
    rxnftrs[wg][nnew].init(1);
    rxnftrs[wg][nnew].set_element(0,e1);
    rxnftrs[wg][nnew].print_atoms();
    nnew++;
  }
#endif
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;

    for (int j=0;j<extent;j++)
    {
      int c1 = climit_l[e1] + j;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      //fflush(stdout);
    }
  } //loop over each element

  nfound = nnew;
//  if (nnew!=nfound) printf(" WARNING: nnew: %i nfound: %i \n",nnew,nfound);

  plus_b = 1;

  printf("\n found %2i features (+b: %i)\n",nfound,plus_b); 

  return nfound;
}


int RXNDB::create_features_climit_ab(int wg)
{
  printf("\n creating element/coordn/add/brk features via climits \n"); 

  int nfound = 0;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;
    nfound += extent;

    printf("  element %2i extent: %i \n",e1,extent);
  }
  nfound *= 2; //for add/brk 
#if NULL_ATOMS
  int nt = NULL_ATOMS;
  nfound += nelemf*nt;
#endif

  if (rxnftrs!=NULL) 
  if (rxnftrs[wg]!=NULL)
  {
    for (int i=0;i<p[wg];i++)
      rxnftrs[wg][i].freemem();
    delete [] rxnftrs[wg];
  }
  rxnftrs[wg] = new RXNFTR[nfound];


 //create feature vectors for each
  int nnew = 0;
#if NULL_ATOMS
  for (int i=0;i<nt*nelemf;i++)
  {
    int e1 = -1;
    rxnftrs[wg][nnew].init(1);
    rxnftrs[wg][nnew].set_element(0,e1);
    rxnftrs[wg][nnew].print_atoms();
    nnew++;
  }
#endif
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;

    for (int j=0;j<extent;j++)
    {
      int c1 = climit_l[e1] + j;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,1,-1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,-1,1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      //fflush(stdout);
    }
  } //loop over each element

  nfound = nnew;
//  if (nnew!=nfound) printf(" WARNING: nnew: %i nfound: %i \n",nnew,nfound);

  plus_b = 1;

  printf("\n found %2i features (+b: %i)\n",nfound,plus_b); 

  return nfound;
}


int RXNDB::create_features_climit_abs(int wg)
{
  printf("\n creating element/coordn/add/brk/swp features via climits \n"); 

  int nfound = 0;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;
    nfound += extent;

    printf("  element %2i extent: %i \n",e1,extent);
  }
  nfound *= 3; //for add/brk/swp
#if NULL_ATOMS
  int nt = NULL_ATOMS;
  nfound += nelemf*nt;
#endif

  if (rxnftrs!=NULL) 
  if (rxnftrs[wg]!=NULL)
  {
    for (int i=0;i<p[wg];i++)
      rxnftrs[wg][i].freemem();
    delete [] rxnftrs[wg];
  }
  rxnftrs[wg] = new RXNFTR[nfound];

  ATOM at1;
  at1.init();
  at1.active = 1;

 //create feature vectors for each
  int nnew = 0;
#if NULL_ATOMS
  for (int i=0;i<nt*nelemf;i++)
  {
    int e1 = -1;
    rxnftrs[wg][nnew].init(1);
    rxnftrs[wg][nnew].set_element(0,e1);
    rxnftrs[wg][nnew].print_atoms();
    nnew++;
  }
#endif
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;

    for (int j=0;j<extent;j++)
    {
      int c1 = climit_l[e1] + j;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,1,-1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,-1,1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,-1,-1,1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      //fflush(stdout);
    }
  } //loop over each element

  nfound = nnew;
//  if (nnew!=nfound) printf(" WARNING: nnew: %i nfound: %i \n",nnew,nfound);

  plus_b = 1;

  printf("\n found %2i features (+b: %i)\n",nfound,plus_b); 

  return nfound;
}

int RXNDB::create_features_climit_abs_att(int wg)
{
  printf("\n creating element/coordn/add/brk/swp features via climits (+attaching/attached) \n"); 

  int nfound = 0;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;
    nfound += extent;

    printf("  element %2i extent: %i \n",e1,extent);
  }
  nfound *= 3; //for add/brk/swp
  nfound += nelemf * nelemf; //attaching pairs
  nfound += nelemf * nelemf; //attached pairs
#if NULL_ATOMS
  int nt = NULL_ATOMS;
  nfound += nelemf*nt;
#endif

  if (rxnftrs!=NULL) 
  if (rxnftrs[wg]!=NULL)
  {
    for (int i=0;i<p[wg];i++)
      rxnftrs[wg][i].freemem();
    delete [] rxnftrs[wg];
  }
  rxnftrs[wg] = new RXNFTR[nfound];

  ATOM at1;
  at1.init();
  at1.active = 1;

 //create feature vectors for each
  int nnew = 0;
#if NULL_ATOMS
  for (int i=0;i<nt*nelemf;i++)
  {
    int e1 = -1;
    rxnftrs[wg][nnew].init(1);
    rxnftrs[wg][nnew].set_element(0,e1);
    rxnftrs[wg][nnew].print_atoms();
    nnew++;
  }
#endif
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;

    for (int j=0;j<nelemf;j++)
    {
      int e2 = emap[j];
      at1.set_element(e2);
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_attaching(0,at1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
    }
    for (int j=0;j<nelemf;j++)
    {
      int e2 = emap[j];
      at1.set_element(e2);
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_attached(0,at1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
    }
    for (int j=0;j<extent;j++)
    {
      int c1 = climit_l[e1] + j;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,1,-1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,-1,1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,-1,-1,1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      //fflush(stdout);
    }
  } //loop over each element

  nfound = nnew;
//  if (nnew!=nfound) printf(" WARNING: nnew: %i nfound: %i \n",nnew,nfound);

  plus_b = 1;

  printf("\n found %2i features (+b: %i)\n",nfound,plus_b); 

  return nfound;
}


int RXNDB::create_features_climit_abs_att_2(int wg)
{
  printf("\n creating element/coordn/add/brk/swp features via climits (+attaching/detaching) \n"); 

  int nex = 0;
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;
    nex += extent;

    printf("  element %2i extent: %i \n",e1,extent);
  }
  int nfound = 3 * nex; //for add/brk/swp
  nfound += nelemf * nelemf; //attaching pairs
  nfound += nelemf * nelemf; //attached pairs
  nfound += nelemf * (nelemf*nelemf); //swp attached/attaching 
#if NULL_ATOMS
  int nt = NULL_ATOMS;
  nfound += nelemf*nt;
#endif

  if (rxnftrs!=NULL) 
  if (rxnftrs[wg]!=NULL)
  {
    for (int i=0;i<p[wg];i++)
      rxnftrs[wg][i].freemem();
    delete [] rxnftrs[wg];
  }
  rxnftrs[wg] = new RXNFTR[nfound];

  ATOM at1,at2;
  at1.init(); at2.init();
  at1.active = 1; at2.active = 1;

 //create feature vectors for each
  int nnew = 0;
#if NULL_ATOMS
  for (int i=0;i<nt*nelemf;i++)
  {
    int e1 = -1;
    rxnftrs[wg][nnew].init(1);
    rxnftrs[wg][nnew].set_element(0,e1);
    rxnftrs[wg][nnew].print_atoms();
    nnew++;
  }
#endif
  for (int i=0;i<nelemf;i++)
  {
    int e1 = emap[i];
    int extent = climit_h[e1] - climit_l[e1] + 1;

    for (int j=0;j<nelemf;j++)
    for (int k=0;k<nelemf;k++)
    {
      int e2 = emap[j];
      int e3 = emap[k];
      at1.set_element(e2);
      at2.set_element(e3);
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_attached(0,at1);
      rxnftrs[wg][nnew].set_attaching(0,at2);
      rxnftrs[wg][nnew].set_abs(0,0,0,1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
    }
    for (int j=0;j<nelemf;j++)
    {
      int e2 = emap[j];
      at1.set_element(e2);
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_attaching(0,at1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
    }
    for (int j=0;j<nelemf;j++)
    {
      int e2 = emap[j];
      at1.set_element(e2);
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_attached(0,at1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
    }
    for (int j=0;j<extent;j++)
    {
      int c1 = climit_l[e1] + j;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,1,-1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,-1,1,-1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      rxnftrs[wg][nnew].init(1);
      rxnftrs[wg][nnew].set_element(0,e1);
      rxnftrs[wg][nnew].set_coordn(0,c1);
      rxnftrs[wg][nnew].set_abs(0,-1,-1,1);
      rxnftrs[wg][nnew].print_atoms();
      nnew++;
      //fflush(stdout);
    }
  } //loop over each element

  nfound = nnew;
//  if (nnew!=nfound) printf(" WARNING: nnew: %i nfound: %i \n",nnew,nfound);

  plus_b = 1;

  printf("\n found %2i features (+b: %i)\n",nfound,plus_b); 

  at1.freemem();
  at2.freemem();

  return nfound;
}

void RXNDB::copy_atom(RXNFTR &ftr1, RXNFTR ftr0)
{
  ATOM at1;
  at1.init();
  at1.active = 1;

  int natoms0 = ftr0.natoms;
  int elem0, coordn0;
  int nadd0,nbrk0,nswp0;
  ftr1.update(natoms0);
  ftr1.reset(natoms0);
  for (int i=0;i<natoms0;i++)
  {
    elem0 = ftr0.get_element(i);
    coordn0 = ftr0.get_coordn(i);
    ftr0.get_abs(i,nadd0,nbrk0,nswp0);
    ftr1.set_element(i,elem0);
    ftr1.set_coordn(i,coordn0);
    ftr1.set_abs(i,nadd0,nbrk0,nswp0);
    if (ftr0.get_nattached(i)>0)
    {
      at1.reset();
      ftr0.get_attached(i,at1);
      ftr1.set_attached(i,at1);
    }
    if (ftr0.get_nattaching(i)>0)
    {
      at1.reset();
      ftr0.get_attaching(i,at1);
      ftr1.set_attaching(i,at1);
    }
  }

  at1.freemem();

  return;
}

void RXNDB::copy_features(int wg1, int wg0)
{
  //printf(" in copy_features: %i %i \n",wg1,wg0);

  RXNFTR* rxnftrs1 = rxnftrs[wg1];
  RXNFTR* rxnftrs0 = rxnftrs[wg0];

  int size1 = p[wg1];
  int size0 = p[wg0];
  copy_features(rxnftrs1,rxnftrs0,size1,size0);

  rxnftrs[wg1] = rxnftrs1; 
  p[wg1] = p[wg0];

  return;
}

void RXNDB::copy_features(RXNFTR* &rxnftrs1, RXNFTR* rxnftrs0, int size1, int size0)
{
  if (rxnftrs1!=NULL)
  { 
    for (int i=0;i<size1;i++)
      rxnftrs1[i].freemem();
    delete [] rxnftrs1;
  }

  rxnftrs1 = new RXNFTR[size0];

  ATOM at1;
  at1.init();
  at1.active = 1;

  for (int i=0;i<size0;i++)
  {
    int natoms0 = rxnftrs0[i].natoms;
    rxnftrs1[i].init(natoms0);
    for (int j=0;j<natoms0;j++)
    {
      int e1 = rxnftrs0[i].get_element(j);
      int c1 = rxnftrs0[i].get_coordn(j);
      int add1,brk1,swp1;
      rxnftrs0[i].get_abs(j,add1,brk1,swp1);
      rxnftrs1[i].set_element(j,e1);
      rxnftrs1[i].set_coordn(j,c1);
      rxnftrs1[i].set_abs(j,add1,brk1,swp1);
      if (rxnftrs0[i].get_nattached(j)>0)
      {
        at1.reset();
        rxnftrs0[i].get_attached(j,at1);
        rxnftrs1[i].set_attached(j,at1);
        //rxnftrs1[i].print_atoms();
      }
      if (rxnftrs0[i].get_nattaching(j)>0)
      {
        at1.reset();
        rxnftrs0[i].get_attaching(j,at1);
        rxnftrs1[i].set_attaching(j,at1);
        //rxnftrs1[i].print_atoms();
      }
    }
    //rxnftrs1[i].print_atoms();
  }
  at1.freemem();

  return;
}

#if 0
int RXNDB::get_index_climit(int ef1, int c1)
{
  int e1 = emap[ef1];
  int index = efindex[ef1] + c1 - climit_l[e1];
  //printf(" get_index, e1: %2i ef1: %i efindex: %i \n",e1,ef1,efindex[ef1]);
  return index; 
}
#endif

void RXNDB::save_structures()
{
  printf("   saving %3i structures \n",N);
  for (int i=0;i<N;i++)
  {
    //write_xyz(i,1,natoms[i],anames[i],xyzr[i],Ea[i]);
    write_xyz(i,-2,natoms[i],anames[i],xyzp[i],Ea[i]);
    write_xyz(i,2,natoms[i],anames[i],xyzp[i],Erxn[i]);
    //write_xyz(i,3,natoms[i],anames[i],xyzt[i],Ea[i]);
  }

  string cmd;
  //cmd = "cat scratch/final.xyzr_* > final.xyzr";
  //system(cmd.c_str());
  cmd = "cat scratch/final.xyz_* > final.xyz";
  system(cmd.c_str());
  cmd = "cat scratch/final.xyzp_* > final.xyzp";
  system(cmd.c_str());
  cmd = "rm scratch/final.xyz*_*";
  system(cmd.c_str());

  return;
}

void RXNDB::write_xyz(int wfile, int type, int natoms1, string* anames1, double* xyz1, double Ea1)
{
  string id = StringTools::int2str(wfile,4,"0");
  string xyzfile_string;
  if (type==1)
    xyzfile_string = "scratch/final.xyzr_"+id;
  else if (type==-2)
    xyzfile_string = "scratch/final.xyz_"+id;
  else if (type==2)
    xyzfile_string = "scratch/final.xyzp_"+id;
  else if (type==3)
    xyzfile_string = "scratch/final.xyzt_"+id;
  else
    return;

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(8);

  xyzfile << natoms1 << endl;
  xyzfile << Ea1 << endl;
  for (int i=0;i<natoms1;i++)
    xyzfile << " " << anames1[i] << " " << xyz1[3*i+0] << " " << xyz1[3*i+1] << " " << xyz1[3*i+2] << endl;

  xyzfile.close();

  return;
}


 //allocs for 1 at a time
void RXNDB::alloc_db(int natoms0)
{
  //if (Nalloc>N+1) return;

  int MAX_ADD = 20;
  int MAX_BRK = 20;

  int N0 = N;
  int N1 = N+1;

  if (N1==1)
  {
    printf("  new database, allocating \n");
    natoms = new int[N1];
    for (int i=0;i<N1;i++)
      natoms[i] = 0;
    anames = new string*[N1];
    for (int i=0;i<N1;i++)
      anames[i] = new string[natoms0];
    anumbers = new int*[N1];
    for (int i=0;i<N1;i++)
      anumbers[i] = new int[natoms0];
    xyzr = new double*[N1];
    for (int i=0;i<N1;i++)
      xyzr[i] = new double[3*natoms0];
    xyzp = new double*[N1];
    for (int i=0;i<N1;i++)
      xyzp[i] = new double[3*natoms0];
    coordnr = new int*[N1];
    for (int i=0;i<N1;i++)
      coordnr[i] = new int[natoms0];
    nadd = new int[N1];
    add = new int*[N1];
    for (int i=0;i<N1;i++)
      add[i] = new int[MAX_ADD];
    nbrks = new int[N1];
    brks = new int*[N1];
    for (int i=0;i<N1;i++)
      brks[i] = new int[MAX_BRK];
    Ea = new double[N1];
    Erxn = new double[N1];

    ids = new int[N1];
    similar = new double[N1];

    qrc = new double*[N1];
    for (int i=0;i<N1;i++)
      qrc[i] = new double[natoms0];
  }
  else
  {
    //printf(" shifting pointers and allocating for %i more structure(s) \n",N1-N0);
    int* natoms1 = new int[N1];
    for (int i=Nalloc;i<N1;i++)
      natoms1[i] = 0;
    for (int i=0;i<N0;i++) natoms1[i] = natoms[i];
    string** anames1 = new string*[N1];
    for (int i=0;i<N0;i++)
      anames1[i] = anames[i];
    for (int i=N0;i<N1;i++)
      anames1[i] = new string[natoms0];
    int** anumbers1 = new int*[N1];
    for (int i=0;i<N0;i++)
      anumbers1[i] = anumbers[i];
    for (int i=N0;i<N1;i++)
      anumbers1[i] = new int[natoms0];
    double** xyzr1 = new double*[N1];
    for (int i=0;i<N0;i++)
      xyzr1[i] = xyzr[i];
    for (int i=N0;i<N1;i++)
      xyzr1[i] = new double[3*natoms0];
    double** xyzp1 = new double*[N1];
    for (int i=0;i<N0;i++)
      xyzp1[i] = xyzp[i];
    for (int i=N0;i<N1;i++)
      xyzp1[i] = new double[3*natoms0];
    int** coordnr1 = new int*[N1];
    for (int i=0;i<N0;i++)
      coordnr1[i] = coordnr[i];
    for (int i=N0;i<N1;i++)
      coordnr1[i] = new int[natoms0];

    int* nadd1 = new int[N1];
    for (int i=0;i<N0;i++) nadd1[i] = nadd[i];
    int** add1 = new int*[N1];
    for (int i=0;i<N0;i++)
      add1[i] = add[i];
    for (int i=N0;i<N1;i++)
      add1[i] = new int[MAX_ADD];

    int* nbrks1 = new int[N1];
    for (int i=0;i<N0;i++) nbrks1[i] = nbrks[i];
    int** brks1 = new int*[N1];
    for (int i=0;i<N0;i++)
      brks1[i] = brks[i];
    for (int i=N0;i<N1;i++)
      brks1[i] = new int[MAX_BRK];

    double* Ea1 = new double[N1];
    for (int i=0;i<N0;i++) Ea1[i] = Ea[i];
    double* Erxn1 = new double[N1];
    for (int i=0;i<N0;i++) Erxn1[i] = Erxn[i];

    int* ids1 = new int[N1];
    for (int i=0;i<N0;i++) ids1[i] = ids[i];
    double* similar1 = new double[N1];
    for (int i=0;i<N0;i++) similar1[i] = similar[i];

    double** qrc1 = new double*[N1];
    for (int i=0;i<N0;i++)
      qrc1[i] = qrc[i];
    for (int i=N0;i<N1;i++)
      qrc1[i] = new double[natoms0];

   //reassign pointers
    delete [] natoms;
    delete [] anames;
    delete [] anumbers;
    delete [] xyzr;
    delete [] xyzp;
    delete [] coordnr;
    delete [] nadd;
    delete [] add;
    delete [] nbrks;
    delete [] brks;
    delete [] Ea;
    delete [] Erxn;
    delete [] ids;
    delete [] similar;
    delete [] qrc;
    natoms = natoms1;
    anames = anames1;
    anumbers = anumbers1;
    xyzr = xyzr1;
    xyzp = xyzp1;
    coordnr = coordnr1;
    nadd = nadd1;
    add = add1;
    nbrks = nbrks1;
    brks = brks1;
    Ea = Ea1;
    Erxn = Erxn1;
    ids = ids1;
    similar = similar1;
    qrc = qrc1;
  }

  Nalloc = N1;

  //printf(" done allocating \n"); fflush(stdout);

  return;
}

void RXNDB::set_limits(int* climit_l0, int* climit_h0)
{
  printf(" saving climits for RXNDB \n");

  climit_l = climit_l0;
  climit_h = climit_h0;

  return;
}

void RXNDB::init(int level0)
{
  quiet = 0;

  level = level0;
  maxg = 100;
  best_fv = new int[maxg];
  for (int i=0;i<maxg;i++) best_fv[i] = 0;
  N = 0;
  p = new int[maxg];
  for (int i=0;i<maxg;i++) p[i] = 0;
  pthresh = new double[maxg];
  for (int i=0;i<maxg;i++) pthresh[i] = PTHRESH0;
  plus_b = 1; //CPMZ just changed
  Nalloc = 0;
  ids = NULL;
  rxnftrs = new RXNFTR*[maxg];
  ga_error = new double[maxg];
  ga_inited = 0;

#if !USE_KNNR
  lsqs = new LSQ[maxg];
  for (int i=0;i<maxg;i++)
    lsqs[i].init();
#else
  knnrs = new KNNR[maxg];
  for (int i=0;i<maxg;i++)
    knnrs[i].init();
#endif

  temperature = 0.;

  printf(" initializing RXNDB at level %i \n",level0);

  natoms = NULL;
  anames = NULL;
  anumbers = NULL;
  xyzr = NULL;
  xyzp = NULL;
  coordnr = NULL;
  nadd = NULL;
  nbrks = NULL;
  add = NULL;
  brks = NULL;
  similar = NULL;
  Ea = NULL;
  Pr = NULL;
  PrT2 = NULL;
  unique = NULL;
  lastAdd = -1;

  nelem = 80;
  emap = NULL;
  remap = NULL;

  naddbrktypes = 0;
  allow_add_brk_elem = NULL;
  addbrk_layer = NULL;
  nrtypes = NULL;
  trxns = NULL;

  nbosr = NULL;
  nbosp = NULL;

  return;
}

void RXNDB::print_reactions()
{
  printf("\n printing reactions from rxn database \n");

  for (int i=0;i<N;i++)
  {
    printf("\n reaction %2i (%i) Ea: %4.1f Erxn: %4.1f \n",i,ids[i],Ea[i],Erxn[i]);

    for (int j=0;j<nadd[i];j++)
      printf(" add: %2i - %2i: %s - %s \n",add[i][2*j+0],add[i][2*j+1],anames[i][add[i][2*j+0]].c_str(),anames[i][add[i][2*j+1]].c_str());
    for (int j=0;j<nbrks[i];j++)
      printf(" brk: %2i x %2i: %s x %s \n",brks[i][2*j+0],brks[i][2*j+1],anames[i][brks[i][2*j+0]].c_str(),anames[i][brks[i][2*j+1]].c_str());
  }

  return;
}


void RXNDB::print_reactions_2()
{
  printf("\n printing reactions from rxn database \n");

  for (int i=0;i<N;i++)
  {
//    printf("\n reaction %2i (%i) Ea: %4.1f Erxn: %4.1f \n",i,ids[i],Ea[i],Erxn[i]);
    printf("\n reaction %4i/%4i Pr: %5.4f nadd: %i nbrk: %i \n",i,ids[i],Pr[i],nadd[i],nbrks[i]);

    printf("  add:");
    for (int j=0;j<nadd[i];j++)
      printf(" %s-%s",anames[i][add[i][2*j+0]].c_str(),anames[i][add[i][2*j+1]].c_str());
    printf("\n");
    printf("  brk:");
    for (int j=0;j<nbrks[i];j++)
      printf(" %s-%s",anames[i][brks[i][2*j+0]].c_str(),anames[i][brks[i][2*j+1]].c_str());
    printf("\n");
  }

  return;
}

int RXNDB::read_xy_data(double* X1, double* Ea1, double* Erxn1, int* ids1, string filename)
{
  printf("\n reading in existing data! \n");

  ifstream output(filename.c_str(),ios::in);
  if (!output)
  {
    printf("   didn't find file %s \n",filename.c_str());
    return 0;
  }
  
  int dim = 40; //hardcoded for now
  string line;
  vector<string> tok_line;

  int nf = 0;
  while (!output.eof())
  {
    for (int i=0;i<5;i++)
    {
      getline(output,line);
      if (output.eof() || (line.size()<1 && i!=4))
      {
        //cout << " RR: " << line << endl;
        printf(" end reached \n");
        break;
      }
      if (i==0)
      {
        //cout << " RR: " << line << endl;
        tok_line = StringTools::tokenize(line, " \t");
        ids1[nf] = atoi(tok_line[3].c_str());
        Ea1[nf] = atof(tok_line[5].c_str());
      //check me
       // Erxn[nf] = atof(tok_line[7].c_str());
      }
      else if (i==3)
      {
        //cout << " RR: " << line << endl;
        tok_line = StringTools::tokenize(line, " \t");
        for (int j=0;j<dim;j++)
          X1[nf*dim+j] = atof(tok_line[j].c_str());
        nf++;
      }
    }
  }

  output.close();

  printf("   found %3i data points \n",nf);

#if 0
  printf("  printing data points \n");
  for (int i=0;i<nf;i++)
  {
    printf("  pt: %4i id: %4i Ea: %6.2f Erxn: %6.2f \n",i,ids1[i],Ea1[i],Erxn1[i]);
    for (int j=0;j<dim;j++)
      printf(" %10.6f",X1[i*dim+j]);
    printf("\n");
  }
#endif

  return nf;
}

void RXNDB::save_xy_data(string filename)
{
  int dim = 40; //hardcoded for now
  double* X = new double[dim*N];
  for (int i=0;i<dim*N;i++) X[i] = 0;
  int* M = new int[dim*N];
  for (int i=0;i<dim*N;i++) M[i] = 0;

#if 0
 //shouldn't need to do this
 //reset any modifications
  for (int i=0;i<N;i++)
  {
    nbosr[i].compare_nbo(nbosp[i],"none",1); //list of broken
    nbosp[i].compare_nbo(nbosr[i],"none",1); //list of added
  }
  fix_unbroken_bonds();
  fix_unadded_bonds();
  reorder_nbo_db();
#endif
  

  int sortme = 0;
  int wg = 0;
  for (int i=0;i<N;i++)
  if (unique[i])
    get_rowX_nbo(wg,i,&X[i*dim],&M[i*dim],sortme);

  ofstream xyfile;
  xyfile.open(filename.c_str());
  xyfile << setprecision(8);

  printf("\n Writing (NBO) X-->Y data to file %s \n",filename.c_str());
  if (sortme==0)
    printf("  this data is NOT sorted or normalized \n");
  else
    printf("  this data is sorted and normalized \n");

  int nu = 0;
  char sbuff[3000];
  for (int i=0;i<N;i++)
  if (unique[i])
  {
    sprintf(sbuff,"  pt: %4i id: %4i Ea: %6.2f Erxn: %6.2f \n",i,ids[i],Ea[i],Erxn[i]); xyfile << sbuff;
    sprintf(sbuff,"   add:"); xyfile << sbuff;
    for (int j=0;j<nadd[i];j++)
    {
      int a1 = add[i][2*j+0];
      int a2 = add[i][2*j+1];
      sprintf(sbuff," %s-%s (%i-%i)",anames[i][a1].c_str(),anames[i][a2].c_str(),a1+1,a2+1); xyfile << sbuff;
    }
    sprintf(sbuff,"\n   brk:"); xyfile << sbuff;
    for (int j=0;j<nbrks[i];j++)
    {
      int b1 = brks[i][2*j+0];
      int b2 = brks[i][2*j+1];
      sprintf(sbuff," %s-%s (%i-%i)",anames[i][b1].c_str(),anames[i][b2].c_str(),b1+1,b2+1); xyfile << sbuff;
    }
    sprintf(sbuff,"\n "); xyfile << sbuff;
    if (M!=NULL)
    {
      for (int j=0;j<dim;j++)
      {
        sprintf(sbuff," %10i",M[i*dim+j]); xyfile << sbuff;
      }
      sprintf(sbuff,"\n "); xyfile << sbuff;
    }
    for (int j=0;j<dim;j++)
    {
      sprintf(sbuff," %10.6f",X[i*dim+j]); xyfile << sbuff;
    }
    //printf("\n");
    xyfile << endl << endl;
    nu++;
  }
  printf("  %4i total data points \n",nu);

  xyfile.close();
  delete [] X;
  if (M!=NULL) delete [] M;

  return;
}


void RXNDB::write_ga(string filename)
{
  printf("  saving to ga file: %s \n",filename.c_str());

  ATOM at1;
  at1.init();

  ofstream file;
  file.open(filename.c_str());
// file << setprecision(8);
// file << fixed;
        
  file << ng << endl << endl;
  
  for (int i=0;i<ng;i++)
  {
    int p1 = p[i];
    int wg = i;

    file << "GA " << i << " features: " << p1 << endl;
    for (int j=0;j<p1;j++)
    {
      int natoms1 = rxnftrs[wg][j].natoms;
      file << " feature " << j << " natoms: " << natoms1 << endl;
      for (int k=0;k<natoms1;k++)
      {
        int wa = k;

        int e1 = rxnftrs[wg][j].get_element(wa);
        int c1 = rxnftrs[wg][j].get_coordn(wa);
        double q = rxnftrs[wg][j].get_q(wa);
        double qsig = rxnftrs[wg][j].get_qsig(wa);
        int a1,b1,s1;
        rxnftrs[wg][j].get_abs(wa,a1,b1,s1);

        int natt = rxnftrs[wg][j].get_nattaching(wa);
        int natd = rxnftrs[wg][j].get_nattached(wa);

        file << "  " << e1 << " " << c1 << " " << q << " " << qsig << " " 
             << a1 << " " << b1 << " " << s1 << " " << endl;
        file << "  " << natt << " " << natd << " " << endl;

        for (int l=0;l<natt;l++)
        {
          rxnftrs[wg][j].get_attaching(wa,at1);
          int e2 = at1.get_element();
          int c2 = at1.get_coordn();
          file << "   " << e2 << " " << c2 << " " << endl;
        }
        for (int l=0;l<natd;l++)
        {
          rxnftrs[wg][j].get_attached(wa,at1);
          int e2 = at1.get_element();
          int c2 = at1.get_coordn();
          file << "   " << e2 << " " << c2 << " " << endl;
        }
      }
    } //loop j over features

    file << endl;
  } //loop i over population

  file.close();

  at1.freemem();

  return;
}

int RXNDB::read_ga(string filename)
{
  //printf(" in read_ga, file: %s \n",filename.c_str()); fflush(stdout);
  int nga = 0;
  int nfound = 0;

  ifstream gafile(filename.c_str(),ios::in);
  if (!gafile)
  {
    printf(" couldn't find GA save file: %s \n",filename.c_str());
    return 0;
  }

  string line;
  vector<string> tok_line;

  getline(gafile,line);
  nga = atoi(line.c_str());
  getline(gafile,line);

  ATOM at1;
  at1.init();
  at1.active = 1;

  int p1 = 0;
  int natt = 0;
  int natd = 0;

  if (nga>0)
  while(!gafile.eof())
  {
    getline(gafile,line);
    //cout << " RR " << line << endl;

    if (line.find("GA")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      p[nfound] = p1 = atoi(tok_line[3].c_str());

      //printf(" p1: %i \n",p1);
      if (p1<1) break;
      rxnftrs[nfound] = new RXNFTR[p1];

      for (int i=0;i<p1;i++)
      {
        getline(gafile,line);
        tok_line = StringTools::tokenize(line, " \t");
        int natoms1 = atoi(tok_line[3].c_str());

        //cout << " RR " << line << endl;
        //printf("  natoms1: %2i \n",natoms1);

        rxnftrs[nfound][i].init(natoms1);
        for (int j=0;j<natoms1;j++)
        {
          int wa = j;
          getline(gafile,line);
          tok_line = StringTools::tokenize(line, " \t");
          int e1 = atoi(tok_line[0].c_str());
          int c1 = atoi(tok_line[1].c_str());
          double q = atof(tok_line[2].c_str());
          double qsig = atof(tok_line[3].c_str());
          int a1 = atoi(tok_line[4].c_str());
          int b1 = atoi(tok_line[5].c_str());
          int s1 = atoi(tok_line[6].c_str());

          getline(gafile,line);
          tok_line = StringTools::tokenize(line, " \t");
          natt = atoi(tok_line[0].c_str());
          natd = atoi(tok_line[1].c_str());

          rxnftrs[nfound][i].set_element(wa,e1);
          rxnftrs[nfound][i].set_coordn(wa,c1);
          rxnftrs[nfound][i].set_abs(wa,a1,b1,s1);
          if (qsig<999.)
          {
            rxnftrs[nfound][i].set_q(wa,q);
            rxnftrs[nfound][i].set_qsig(wa,qsig);
          }

          int e2 = -1;
          int c2 = -1;
          for (int k=0;k<natt;k++)
          {
            getline(gafile,line);
            tok_line = StringTools::tokenize(line, " \t");
            e2 = atoi(tok_line[0].c_str());
            c2 = atoi(tok_line[1].c_str());
 
            at1.set_element(e2);
            at1.set_coordn(c2);

            rxnftrs[nfound][i].set_attaching(wa,at1);
          }
          for (int k=0;k<natd;k++)
          {
            getline(gafile,line);
            tok_line = StringTools::tokenize(line, " \t");
            e2 = atoi(tok_line[0].c_str());
            c2 = atoi(tok_line[1].c_str());
 
            at1.set_element(e2);
            at1.set_coordn(c2);

            rxnftrs[nfound][i].set_attached(wa,at1);
          }
        } //loop j over atoms in features
      } //loop i over features

      nfound++;
    } //if found GA section

  } //if !eof

  gafile.close();

  if (nga!=nfound)
  {
    printf("  problem reading %s, not enough data \n",filename.c_str());
    exit(1);
  }

#if 0
  for (int i=0;i<nfound;i++)
  {
    printf("\n features from file for GA %2i \n",i);
    for (int j=0;j<p[i];j++)
      rxnftrs[i][j].print_atoms();
  }
  printf("\n");
#endif

  at1.freemem();

  return nfound;
}

int RXNDB::find_ts(int ws)
{
#if !TEST_GA
  return -1;
#endif

  int wsdb = -1;
  for (int i=0;i<N;i++)
  if (ids[i]==ws)
  {
    wsdb = i;
    break;
  }
  return wsdb;
}


int RXNDB::get_abtype(int i, int nadd1, int nbrks1)
{
  int abtype = get_abtype(nadd1,nbrks1);
  if (abtype==-1)
  {
    printf("\n bad abtype (%i). id: %4i ids: %4i \n",abtype,i,ids[i]);
    printf(" exiting! \n");
    exit(1);
  }

  return abtype;
}

int RXNDB::get_abtype(int nadd1, int nbrks1)
{
  //could change to square

  int abtype = -1;
  if (nadd1==1 && nbrks1==0)
    abtype = 0;
  else if (nadd1==0 && nbrks1==1)
    abtype = 1;
  else if (nadd1==2 && nbrks1==0)
    abtype = 2;
  else if (nadd1==0 && nbrks1==2)
    abtype = 3;
  else if (nadd1==1 && nbrks1==1)
    abtype = 4;
  else if (nadd1==2 && nbrks1==1)
    abtype = 5;
  else if (nadd1==1 && nbrks1==2)
    abtype = 6;
  else if (nadd1==2 && nbrks1==2)
    abtype = 7;
  else if (nadd1==3 && nbrks1==2)
    abtype = 8;
  else if (nadd1==2 && nbrks1==3)
    abtype = 9;
  else if (nadd1==3 && nbrks1==3)
    abtype = 10;
  else if (nadd1>=3 || nbrks1>=3)
    abtype = 11;


  return abtype;
}

void RXNDB::set_unique(int N1, int* unique1, int* ids1)
{
  if (N1<N)
  {
    printf("\n  ERROR: unique data points array wrong size: %i vs %i \n",N1,N);
    exit(1);
  }

  if (unique!=NULL) delete [] unique;
  unique = new int[N];

  int wid = 0;
  for (int i=0;i<N;i++)
  for (int j=0;j<N1;j++)
  if (ids[i]==ids1[j])
  {
    unique[i] = unique1[j];
    wid++;
  }

  printf("\n  added %2i unique data point labels \n",wid);
  if (wid!=N)
  {
    printf(" ERROR: unique data label size not correct. wid: %2i N: %2i \n",wid,N);
    exit(1);
  }


  return;
}


void RXNDB::add_nbo_data(int nnbo, NBO* nbo1r, NBO* nbo1p, int* ids1)
{
  if (ids==NULL)
  { 
    printf(" cannot add NBO data, ids problem \n");
    return;
  }
  if (ids1==NULL)
  { 
    printf(" cannot add NBO data, ids1 problem \n");
    return;
  }
  if (N<1)
  {
    printf(" RXNDB not ready, cannot add NBO data \n");
    return;
  }
  if (nbosr!=NULL) delete [] nbosr;
  nbosr = new NBO[N];
  if (nbosp!=NULL) delete [] nbosp;
  nbosp = new NBO[N];

  printf(" adding NBO data to RXNDB. nnbo: %i N: %i \n",nnbo,N);
  //for (int i=0;i<min(nnbo,N);i++) printf(" ids: %4i %4i \n",ids[i],ids1[i]);

  int wid = 0;
  for (int i=0;i<N;i++)
  for (int j=0;j<nnbo;j++)
  if (ids[i]==ids1[j])
  {
    nbosr[i] = nbo1r[j];
    nbosp[i] = nbo1p[j];
    wid++;
  }

  printf("  added %2i NBO data points \n",wid);
  if (wid!=N)
  {
    printf(" ERROR: NBO sizes not correct. wid: %2i N: %2i \n",wid,N);
    exit(1);
  }

  fix_unbroken_bonds();
  fix_unadded_bonds();
  reorder_nbo_db();

  return;
}

void RXNDB::fix_unbroken_bonds()
{
  ICoord ic1;
  ic1.alloc(500);

 //check that breaks are in broken bond list of nbosr
  //adds broken to end of list
  for (int i=0;i<N;i++)
  {
    for (int j=0;j<nbrks[i];j++)
    {
      int bb1 = brks[i][2*j+0];
      int bb2 = brks[i][2*j+1];

      //printf(" bb1/2: %2i %2i",bb1+1,bb2+1);
      int found = 0;
      for (int k=0;k<nbosr[i].nb;k++)
      {
        int b1  = nbosr[i].blist[k];
        int bt1 = nbosr[i].bmo_atoms[2*b1+0]-1;
        int bt2 = nbosr[i].bmo_atoms[2*b1+1]-1;
        //printf(" bt: %2i %2i",bt1+1,bt2+1);
        if (bt1==bb1 && bt2==bb2)
          found = 1;
        if (bt1==bb2 && bt2==bb1)
          found = 1;
        if (found) break;
      }
      //printf("\n");
      if (!found)
      {
        ic1.reset(natoms[i],anames[i],anumbers[i],xyzp[i]);
       //instead of one half this distance??
        double d0 = (ic1.getR(bb1) + ic1.getR(bb2))/2.0; 
        double d1 = ic1.distance(bb1,bb2);
        if (d1>d0)
        {
          printf("    id: %4i connection broken: %2i-%2i but no bond broken. distance: %5.2f/%5.2f \n",ids[i],bb1+1,bb2+1,d1,d0);
          //nbosr[i].print_nbo();
          //nbosr[i].compare_nbo(nbosp[i],"none",0);
          int bnew = -1;
          for (int k=0;k<nbosr[i].bmo;k++)
          {
            int bn1 = nbosr[i].bmo_atoms[2*k+0]-1;
            int bn2 = nbosr[i].bmo_atoms[2*k+1]-1;
            if (bn1==bb1 && bn2==bb2)
              bnew = k;
            if (bn1==bb2 && bn2==bb1)
              bnew = k;
            if (bnew!=-1) break;
          }
          if (bnew>-1)
          {
            printf("     adding %2i to blist (nb: %i) \n",bnew,nbosr[i].nb);
            nbosr[i].blist[nbosr[i].nb] = bnew;
            nbosr[i].nb++;
          }
        }
      }
    } //loop j over nbrks
  } //loop i over N

  ic1.freemem();

  return;
}

void RXNDB::fix_unadded_bonds()
{
  ICoord ic1;
  ic1.alloc(500);

 //check that adds in nbosp are real
  for (int i=0;i<N;i++)
  {
    ic1.reset(natoms[i],anames[i],anumbers[i],xyzp[i]);
    for (int j=0;j<nbosp[i].nb;j++)
    {
      int a1  = nbosp[i].blist[j];
      int at1 = nbosp[i].bmo_atoms[2*a1+0]-1;
      int at2 = nbosp[i].bmo_atoms[2*a1+1]-1;
      if (at2<0) at2 = at1;

      if (at1<0 || at2<0) { printf(" ERROR on ids: %4i \n",ids[i]); exit(1); }
     //instead of one half this distance??
      double d0 = (ic1.getR(at1) + ic1.getR(at2))/2.0;
      double d1 = 0.;
      if (at1!=at2) 
        d1 = ic1.distance(at1,at2);
      if (d1>d0)
      {
        printf("    id: %4i nbo added, but not connected: %2i-%2i. distance: %5.2f/%5.2f \n",ids[i],at1+1,at2+1,d1,d0);
        for (int k=j;k<nbosp[i].nb-1;k++)
          nbosp[i].blist[k] = nbosp[i].blist[k+1];
        nbosp[i].nb--;
      }
    } //loop j over nbosp[i].nb (added orbs)
  } //loop i over N

  ic1.freemem();

  return;
}

void RXNDB::reorder_nbo_db()
{
  printf("   reordering blist to match add/brk order \n"); 

  for (int i=0;i<N;i++)
  {
    int shift = 0;
    for (int j=0;j<nbosp[i].nb - shift;j++)
    {
      int a1  = nbosp[i].blist[j];
      int at1 = nbosp[i].bmo_atoms[2*a1+0]-1;
      int at2 = nbosp[i].bmo_atoms[2*a1+1]-1;

      int found = 0;
      for (int k=0;k<nadd[i];k++)
      {
        int aa1 = add[i][2*k+0];
        int aa2 = add[i][2*k+1];
        if (at1==aa1 && at2==aa2)
          found = 1;
        if (at1==aa2 && at2==aa1)
          found = -1;

//CPMZ check me
        if (found==-1)
        {
          nbosp[i].bmo_atoms[2*a1+0] = at2+1;
          nbosp[i].bmo_atoms[2*a1+1] = at1+1;
        }
      }
      if (!found)
      {
        printf("  rxn %4i id: %4i  add %2i-%2i not in connectivity list \n",i,ids[i],at1+1,at2+1);
        int t1 = nbosp[i].nb-1;
        int tt1 = nbosp[i].blist[j];
        for (int k=j;k<t1;k++)
          nbosp[i].blist[k] = nbosp[i].blist[k+1];
        nbosp[i].blist[t1] = tt1;
        shift++;
        j--;
      }
    } //loop j over nb
  } //loop i over N

  for (int i=0;i<N;i++)
  {
    int shift = 0;
    for (int j=0;j<nbosr[i].nb - shift;j++)
    {
      int b1  = nbosr[i].blist[j];
      int bt1 = nbosr[i].bmo_atoms[2*b1+0]-1;
      int bt2 = nbosr[i].bmo_atoms[2*b1+1]-1;

      int found = 0;
      for (int k=0;k<nbrks[i];k++)
      {
        int bb1 = brks[i][2*k+0];
        int bb2 = brks[i][2*k+1];
        if (bt1==bb1 && bt2==bb2)
          found = 1;
        if (bt1==bb2 && bt2==bb1)
          found = -1;

        if (found==-1)
        {
          nbosr[i].bmo_atoms[2*b1+0] = bt2+1;
          nbosr[i].bmo_atoms[2*b1+1] = bt1+1;
        }
      }
      if (!found)
      {
        printf("  rxn %4i id: %4i  brk %2i-%2i not in connectivity list \n",i,ids[i],bt1+1,bt2+1);
        int t1 = nbosr[i].nb-1;
        int tt1 = nbosr[i].blist[j];
        for (int k=j;k<t1;k++)
          nbosr[i].blist[k] = nbosr[i].blist[k+1];
        nbosr[i].blist[t1] = tt1;
        shift++;
        j--;
      }
    } //loop j over nb
  } //loop i over N

  return;
}

void RXNDB::find_pthresh(int wg, double pfalseneg)
{
  

  return;
}


