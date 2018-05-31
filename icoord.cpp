// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "icoord.h"
#include "utils.h"
using namespace std;

#define MAX_FRAG_DIST 12.0


//NOTE: eta4,5,6 --> break into eta2's


//supports up to 1 TM 
int ICoord::ic_create_tm()
{
  printf("  in ic_create_tm \n");

 //note: split_h called from ic_create makes H one coordinate 
 //       unless TM or B bridging
 // print_bonds();

  int tbond = 0;
  int tbond2 = 0;

  nbondstm = 0;
  for (int i=0;i<25;i++)
  for (int j=0;j<15;j++)
    bondstm[i][j] = -1;
  nbondstm2 = 0;
  for (int i=0;i<25;i++)
  for (int j=0;j<15;j++)
    bondstm2[i][j] = -1;

  int wtm = -1;
 //get everything bonded to first TM
  for (int i=0;i<natoms;i++)
  if (isTM(i))
  {
    printf("   found TM: %i/%s \n",i+1,anames[i].c_str());
    wtm = i;
    for (int j=0;j<nbonds;j++)
    {
      if (bonds[j][0]==wtm)
      {
        bondstm[tbond][0] = wtm;
        bondstm[tbond++][1] = bonds[j][1];
      }
      else if (bonds[j][1]==wtm)
      {
        bondstm[tbond][0] = wtm;
        bondstm[tbond++][1] = bonds[j][0];
      }
    }
    break;
  } //loop i over TM atoms

  int wtm2 = -1;
 //get everything bonded to second TM
  for (int i=0;i<natoms;i++)
  if (isTM(i) && wtm!=i)
  {
    printf("   found TM: %i/%s \n",i+1,anames[i].c_str());
    wtm2 = i;
    for (int j=0;j<nbonds;j++)
    {
      if (bonds[j][0]==wtm2)
      {
        bondstm2[tbond2][0] = wtm2;
        bondstm2[tbond2++][1] = bonds[j][1];
      }
      else if (bonds[j][1]==wtm2)
      {
        bondstm2[tbond2][0] = wtm2;
        bondstm2[tbond2++][1] = bonds[j][0];
      }
    }
    break;
  } //loop i over TM atoms

  if (wtm2<0)
    printf("   no second TM here \n");

#if 1
  printf("  printing raw TM bonds for first TM \n");
  for (int i=0;i<tbond;i++)
  {
    printf("  ");
    for (int j=0;j<15;j++)
    if (bondstm[i][j]>-1)
      printf(" %2i ",bondstm[i][j]+1);
    printf("\n");
  }
  if (tbond2>0)
    printf("  printing raw TM bonds for second TM \n");
  for (int i=0;i<tbond2;i++)
  {
    printf("  ");
    for (int j=0;j<15;j++)
    if (bondstm2[i][j]>-1)
      printf(" %2i ",bondstm2[i][j]+1);
    printf("\n");
  }
#endif

  get_eta(tbond,bondstm);
  get_eta(tbond2,bondstm2);

 //count total # attached
  coordntm = 0;
  for (int i=0;i<tbond;i++)
  for (int j=1;j<15;j++)
  if (bondstm[i][j]>-1)
    coordntm++;
  coordntm2 = 0;
  for (int i=0;i<tbond2;i++)
  for (int j=1;j<15;j++)
  if (bondstm2[i][j]>-1)
    coordntm2++;

  nbondstm = tbond;
  nbondstm2 = tbond2;

  printf("  coordntm: %2i nbondstm: %2i \n",coordntm,nbondstm);
  if (wtm2>-1)
    printf("  coordntm2: %2i nbondstm2: %2i \n",coordntm2,nbondstm2);

  return tbond;
}


void ICoord::get_eta(int& tbond, int** btm)
{
  //printf("  in get_eta for tbond: %i \n",tbond);

 //collect eta-2+ ligands
  for (int i=0;i<tbond;i++)
  for (int j=1;j<15;j++)
  if (btm[i][j]>-1)
  {
    for (int k=0;k<i;k++)
    for (int l=1;l<15;l++)
    if (btm[k][l]>-1)
    {

      int a1 = btm[i][j];
      int a2 = btm[k][l];
      //printf("  i,k: %i %i a1,a2: %2i %2i \n",i,k,a1+1,a2+1);
      if (bond_exists(a1,a2))
      {
        printf("   atoms eta to metal: %2i %2i \n",a1+1,a2+1);
        int nf = 0;
        int f2 = 0;
        for (int p=1;p<15;p++)
        if (btm[i][j+nf]<0)
        {
          //printf("    p: %i btm[k][p]: %i \n",p,btm[k][p]+1);
          f2 = 0;
          for (int q=1;q<15;q++)
          if (btm[k][p]==btm[i][q])
            f2 = 1;
          if (!f2)
          {
            //printf("     %i not found, adding \n",btm[k][p]+1);
            btm[i][j+nf++] = btm[k][p];
          }
        }
        else
        {
          nf++;
          p--;
        }
        if (f2>0)
        for (int p=0;p<15;p++)
          btm[k][p] = -1;
      } //if bond exists

    } //loop k < i
  } //loop i over nbtm

 //rearrange btm so no entries are void
  for (int i=0;i<tbond;i++)
  {
    int f1 = 0;
    for (int j=0;j<15;j++)
    if (btm[i][j]>-1)
      f1++;
    if (f1==0)
    {
      for (int j=i;j<tbond-1;j++)
      for (int k=0;k<15;k++)
        btm[j][k] = btm[j+1][k];
      tbond--;
      i--; //in case next one is void also
    }
  }

#if 1
  printf("   printing reorg TM bonds \n");
  for (int i=0;i<tbond;i++)
  {
    printf("  ");
    for (int j=0;j<15;j++)
    if (btm[i][j]>-1)
      printf(" %2i ",btm[i][j]+1);
    printf("\n");
  }
#endif

  for (int i=0;i<tbond;i++)
  {
    int ntb = 0;
    for (int j=0;j<15;j++)
    if (btm[i][j]>-1)
      ntb++;

    //printf("   tbond %i connections: %i \n",i,ntb-1);
    if (ntb==3) //eta2 type
    {
      //currently: keep "eta2" for bond, counts as 1 coord
      int tm1 = btm[i][0];
      coordn[tm1]--;
    }
    else if (ntb==4) //eta3 type
    {
     //deletes central connection in eta3, counts as 2 coord
      int tm1 = btm[i][0];
      coordn[tm1]--;

      int b1  = btm[i][1];
      int b2  = btm[i][2];
      int b3  = btm[i][3];
      int wb = -1; int ib = -1;
      if (bond_exists(b1,b2) && bond_exists(b2,b3)) { wb = b2; ib = 2; }
      if (bond_exists(b1,b3) && bond_exists(b2,b3)) { wb = b3; ib = 3; }
      if (bond_exists(b1,b2) && bond_exists(b1,b3)) { wb = b1; ib = 1; }
      printf("   breaking middle link: tm-%i \n",wb+1);
      if (ib==1)
      {
        btm[i][1] = b2;
        btm[tbond][1] = b3;
      }
      else if (ib==2)
      {
        btm[i][1] = b1;
        btm[tbond][1] = b3;
      }
      else if (ib==3)
      {
        btm[i][1] = b1;
        btm[tbond][1] = b2;
      }
      btm[tbond][0] = tm1;
      for (int j=2;j<15;j++) btm[i][j] = -1;
      for (int j=2;j<15;j++) btm[tbond][j] = -1;
      tbond++;
      coordn[wb]--;

     //old, just keep 2 links
//      for (int j=ib;j<3;j++)
//        btm[i][j] = btm[i][j+1];
//      btm[i][3] = -1;
    }
    else if (ntb==5)
    {
      printf("\n   TESTING eta4 \n");
      int tm1 = btm[i][0];
      coordn[tm1] -= 2;

     //find most distant pair
      double dmax = 0.;
      int wb1 = -1; int wb2 = -1;
      for (int j=1;j<5;j++)
      {
        int a1 = btm[i][j];
        for (int k=1;k<j;k++)
        { 
          int a2 = btm[i][k];
          double d1 = distance(a1,a2);
          if (d1>dmax) 
          {
            dmax = d1;
            wb1 = a1;
            wb2 = a2;
          }
        }
      }
      printf("   distant pair: %i %i \n",wb1+1,wb2+1);

      btm[i][1] = wb1;
      btm[tbond][1] = wb2;
      btm[tbond][0] = tm1;
      for (int j=2;j<15;j++) btm[i][j] = -1;
      for (int j=2;j<15;j++) btm[tbond][j] = -1;
      tbond++;
    }
    else if (ntb>5)
    {
      printf("\n\n ERROR: eta > 4 not supported \n");
//      exit(1);
    }
  }

//  printf(" ENDING EARLY \n"); exit(1);

  for (int i=0;i<tbond;i++)
  {
    int tb1 = 0;
    for (int j=1;j<15;j++)
    if (btm[i][j]>-1)
      tb1++;

    for (int j=0;j<i;j++)
    {
      int tb2 = 0;
      for (int k=1;k<15;k++)
      if (btm[j][k]>-1)
        tb2++;
      int tswp = 0;
      if (tb1>tb2)
        tswp = 1;
      if (tb1==tb2 && btm[i][tb1]<btm[j][tb2])
        tswp = 1;
      if (tswp)
      {
        int* tbtmp = btm[i];
        btm[i] = btm[j];
        btm[j] = tbtmp;
        i--; 
        break;
      }
    }
  }

#if 1
  printf("   printing reconfig TM bonds \n");
  for (int i=0;i<tbond;i++)
  {
    printf("  ");
    for (int j=0;j<15;j++)
    if (btm[i][j]>-1)
      printf("  %2i",btm[i][j]+1);
    printf("\n");
  }
#endif

  return;
}


int ICoord::init(string xyzfile){


// printf(" xyzfile: %s \n",xyzfile);
 //printf("\n");
 //cout << " xyzfile: " << xyzfile << endl;
 structure_read(xyzfile);
 
 //print_xyz();

 alloc_mem();
 // printf(" done allocating memory\n");

 farBond = 1.0;
 int done = ic_create();

// printf(" initializing MM parameters \n");
 mm_init();

 //printf("\n\n");

 return 1;
}



// initialize by feeding in xyz coordinates
int ICoord::init(int nat, string* anam, int* anum, double* xyz){

// printf(" initializing icoord via xyz structure \n");
 natoms = nat;
// printf(" natoms: %i \n",nat);
// for (int i=0;i<natoms;i++)
//    printf(" %s %1.3f %1.3f %1.3f \n",anam[i].c_str(),xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]);
// for (int i=0;i<natoms;i++)
//    printf(" %2i\n",anum[i]);
// fflush(stdout);

//otherwise allocated in structure_read
 anumbers = new int[1+natoms];
 amasses = new double[1+natoms];
 anames = new string[1+natoms];
 coords = new double[natoms*3];
 //coordsts = new double[natoms*3];
 coords0 = new double[natoms*3];

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

 alloc_mem();

 farBond = 1.0;
 int done = ic_create();


// printf(" initializing MM parameters \n");
 mm_init();

 //printf("\n\n");

 return 1;
}

// initialize memory only
int ICoord::alloc(int size){

 natoms = size;

//otherwise allocated in structure_read
 anumbers = new int[1+natoms];
 amasses = new double[1+natoms];
 anames = new string[1+natoms];
 coords = new double[natoms*3];
 //coordsts = new double[natoms*3];
 coords0 = new double[natoms*3];

 alloc_mem();

 farBond = 1.0;

 return 1;
}


// initialize by feeding in xyz coordinates
int ICoord::reset(double* xyz){

// printf(" resetting icoord via xyz structure \n");

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

// int done = ic_create();

// printf(" initializing MM parameters \n");
// mm_init();

 //printf("\n\n");

 return 1;
}

// initialize by feeding in xyz coordinates
int ICoord::reset(int nat, string* anam, int* anum, double* xyz){

// printf(" resetting icoord via xyz structure \n");
 natoms = nat;

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

// int done = ic_create();

// printf(" initializing MM parameters \n");
// mm_init();

 //printf("\n\n");

 return 1;
}

void ICoord::update_ic(){

  update_bonds();
  update_angles();
  update_torsion();
  update_imptor();
  update_nonbond();

  return;
} 
 
void ICoord::create_xyz()
{
  double* nxyz = new double[3*natoms];
  printf ("xyz_create not implemented\n");
  int* adone = new int[natoms];
  for (int i=0;i<natoms;i++) adone[i]=0;
  adone[0]=1;
  nxyz[0] = coords[0];
  nxyz[1] = coords[1];
  nxyz[2] = coords[2];
  
  double* v1 = new double[3];
  double* u1 = new double[3];
  v1[0] = coords[3] - coords[0];
  v1[1] = coords[4] - coords[1];
  v1[2] = coords[5] - coords[2];

  double norm1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  u1[0]=v1[0]/norm1;
  u1[1]=v1[1]/norm1;
  u1[2]=v1[2]/norm1;

  double R;
  for (int i=0;i<nbonds;i++)
    if ((bonds[i][0]==0 && bonds[i][1]==1) ||
        (bonds[i][1]==0 && bonds[i][0]==1) )
      R = bondd[i];

//alternatively, just put this along the x axis
  nxyz[3] = u1[0]*R;
  nxyz[4] = u1[1]*R;
  nxyz[5] = u1[2]*R;

  for (int i=0;i<nangles;i++)
  {

  }
 
  delete [] adone;

  return;
}


int ICoord::ic_create()
{
  //printf(" Creating internals from xyz \n"); fflush(stdout);
  make_bonds();
  split_h();
  coord_num(); // counts # surrounding species

  //printf(" now making angles \n"); fflush(stdout);
  make_angles();
  //printf(" now making torsions \n"); fflush(stdout);
  make_torsions();
  //printf(" now making improper torsions \n"); fflush(stdout);
  make_imptor();

  //printf(" now counting nonbond\n"); fflush(stdout);
  n_nonbond = make_nonbond(); //anything not connected by bond or angle

//  print_ic();

  return 0;
}

int ICoord::ic_create_nobonds()
{
//  printf(" Creating internals from xyz, skipping bond making \n");
  coord_num(); // counts # surrounding species
  make_angles();
  make_torsions();
  make_imptor_nobonds();
  n_nonbond = make_nonbond(); //anything not connected by bond or angle

//  print_ic();

  return 0;
}


void ICoord::make_frags() 
{
  //printf(" in make_frags() \n");

  frags = new int[natoms];
  nfrags = 0;

  for (int i=0;i<natoms;i++) frags[i] = -1;

  for (int i=0;i<nbonds;i++)
  {
    if ( frags[bonds[i][0]] == -1 && frags[bonds[i][1]] == -1 )
    {
      frags[bonds[i][0]] = nfrags;
      frags[bonds[i][1]] = nfrags;
      nfrags++;
    }
    else if ( frags[bonds[i][0]] == -1 && frags[bonds[i][1]] > -1 )
    {
      frags[bonds[i][0]] = frags[bonds[i][1]];
    }
    else if ( frags[bonds[i][0]] > -1 && frags[bonds[i][1]] == -1 )
    {
      frags[bonds[i][1]] = frags[bonds[i][0]];
    }
    else if ( frags[bonds[i][0]] > -1 && frags[bonds[i][1]] > -1 )
    {
      //both frags assigned already
      if (frags[bonds[i][0]] != frags[bonds[i][1]])
      {
        //printf(" WARNING, need to merge: %i %i \n",frags[bonds[i][0]],frags[bonds[i][1]]);
        //mergelist[2*nmerge+0] = frags[bonds[i][0]];
        //mergelist[2*nmerge+1] = frags[bonds[i][1]];
        int f1 = min(frags[bonds[i][0]],frags[bonds[i][1]]);
        int f2 = max(frags[bonds[i][0]],frags[bonds[i][1]]);
        //printf(" merging: %i/%i \n",f1,f2);
        for (int j=0;j<natoms;j++)
        if (frags[j]==f2)
          frags[j] = f1;
        if (f2==nfrags-1)
          nfrags--;
      }
    }
  } //loop i over nbonds

  for (int i=0;i<natoms;i++)
  if (frags[i]==-1)
    frags[i] = nfrags++;

//  for (int i=0;i<natoms;i++)
//    printf(" atom[%i] frag: %i \n",i,frags[i]);

  return;
}

void ICoord::bond_frags() 
{
//  printf(" in bond_frags() \n");
  if (nfrags<2) return;

  int found = 0;
  int found2 = 0;
  int found3 = 0;

  int a1,a2;
  int b1,b2;
  int c1,c2;
  double mclose;
  double mclose2;
  double mclose3;
  for (int n1=0;n1<nfrags;n1++)
  for (int n2=0;n2<n1;n2++)
  {
//    if (natoms<50)
//      printf(" connecting frag %i to frag %i: ",n1,n2);

    found = 0;
    found2 = 0;
    found3 = 0;
    double close = 0.;
    mclose = 1000.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<natoms;j++)
    if (frags[i]==n1 && frags[j]==n2)
    {
      close = distance(i,j);
      if (close<mclose && close < MAX_FRAG_DIST)
      {
        mclose = close;
        a1 = i;
        a2 = j;
        found = 1;
      }
    }

#if 0   
   //connect second pair, heavies or H-bond only, away from 1st pair
    b1 = -1;
    b2 = -1;
    mclose2 = 1000.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<natoms;j++)
    if (frags[i]==n1 && frags[j]==n2)
    {
      close = distance(i,j);
      double dia1 = distance(i,a1);
      double dja1 = distance(j,a1);
      double dia2 = distance(i,a2);
      double dja2 = distance(j,a2);
      double dist21 = (dia1+dja1)/2.;
      double dist22 = (dia2+dja2)/2.;
      if (anumbers[i] > 1 || anumbers[j] > 1)
      if (dist21 > 4.5 && dist22 > 4.5) //standard
//      if (dia1 > 4.5 && dja1 > 4.5 && dia2 > 4.5 && dja2 > 4.5) //possible change
      if (close<mclose2 && close < MAX_FRAG_DIST)
      {
        mclose2 = close;
        b1 = i;
        b2 = j;
        found2 = 1;
      }
    }

   //connect third pair, heavies or H-bond only, away from 1st pair
    c1 = -1;
    c2 = -1;
    mclose3 = 1000.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<natoms;j++)
    if (frags[i]==n1 && frags[j]==n2)
    {
      close = distance(i,j);
      double dia1 = distance(i,a1);
      double dja1 = distance(j,a1);
      double dia2 = distance(i,a2);
      double dja2 = distance(j,a2);
      double dib1 = distance(i,b1);
      double djb1 = distance(j,b1);
      double dib2 = distance(i,b2);
      double djb2 = distance(j,b2);
      double dist31 = (dia1+dja1)/2.;
      double dist32 = (dia2+dja2)/2.;
      double dist33 = (dib1+djb1)/2.;
      double dist34 = (dib2+djb2)/2.;
      if (anumbers[i] > 1 || anumbers[j] > 1)
      if (dist31 > 4.5 && dist32 > 4.5 && dist33 > 4.5 && dist34 > 4.5) //standard
//      if (dia1 > 4.5 && dja1 > 4.5 && dia2 > 4.5 && dja2 > 4.5)
//      if (dib1 > 4.5 && djb1 > 4.5 && dib2 > 4.5 && djb2 > 4.5)
      if (close<mclose3 && close < MAX_FRAG_DIST)
      {
        mclose3 = close;
        c1 = i;
        c2 = j;
        found3 = 1;
      }
    }
#endif

    if (found && !bond_exists(a1,a2))
    {
      //printf(" bond pair1 added : %i %i ",a1,a2);
      bonds[nbonds][0] = a1;
      bonds[nbonds][1] = a2;
      bondd[nbonds] = mclose;
      nbonds++;
    } // if found

#if 0
    if (found2 && !bond_exists(b1,b2))
    {
      printf(" bond pair2 added : %i %i ",b1,b2);
      bonds[nbonds][0] = b1;
      bonds[nbonds][1] = b2;
      bondd[nbonds] = mclose2;
      nbonds++;
    } // if found

    if (found3 && !bond_exists(c1,c2))
    {
      printf(" bond pair3 added : %i %i ",c1,c2);
      bonds[nbonds][0] = c1;
      bonds[nbonds][1] = c2;
      bondd[nbonds] = mclose3;
      nbonds++;
    } // if found
#endif

//    if (found || natoms<50)
//      printf("\n");
  }//loop n over nfrags



  return;
}

void ICoord::update_bonds(){  
  for (int i=0;i<nbonds;i++)
    bondd[i] = distance(bonds[i][0],bonds[i][1]);
  return;
}

void ICoord::update_angles(){
  for (int i=0;i<nangles;i++)
    anglev[i] = angle_val(angles[i][0],angles[i][1],angles[i][2]);
  return;
}

void ICoord::update_torsion(){
  for (int i=0;i<ntor;i++)
    torv[i]=torsion_val(torsions[i][0],torsions[i][1],torsions[i][2],torsions[i][3]);
  return;
}

void ICoord::update_imptor(){
  for (int i=0;i<nimptor;i++)
    imptorv[i]=torsion_val(imptor[i][0],imptor[i][1],imptor[i][2],imptor[i][3]);
  return;
}

void ICoord::update_nonbond(){
  for (int i=0;i<n_nonbond;i++)
    nonbondd[i] = distance(nonbond[i][0],nonbond[i][1]);
  return;
}

void ICoord::split_h()
{
  coord_num();

  for (int i=0;i<natoms;i++)
  if (anumbers[i]==1 && coordn[i]>1)
  { 
    //printf(" found %i coordinate H: %i \n",coordn[i],i+1);

    for (int j=0;j<nbonds;j++)
    if (bonds[j][0]==i || bonds[j][1]==i)
    {
      //printf("  connected to %i \n",j+1);
      int b1 = -1; //bond to delete
      for (int k=0;k<j;k++)
      if (b1==-1)
      if (bonds[k][0]==i || bonds[k][1]==i)
      {
        b1 = j;
        if (bondd[k]>bondd[j])
        {
          //printf(" bond %i is longer than %i (1) \n",k+1,j+1);
          b1 = k;
        }
        //else
          //printf(" bond %i is longer than %i (2) \n",j+1,k+1);
        if (anumbers[bonds[j][0]]==5 || anumbers[bonds[j][1]]==5
         || anumbers[bonds[k][0]]==5 || anumbers[bonds[k][1]]==5)
          b1 = -2;
        if (isTM(bonds[j][0]) || isTM(bonds[j][1])
         || isTM(bonds[k][0]) || isTM(bonds[k][1]))
          b1 = -2;

        if (b1>-1)
        {
          for (int l=b1;l<nbonds-1;l++)
          {
            bonds[l][0] = bonds[l+1][0];
            bonds[l][1] = bonds[l+1][1];
            bondd[l]    = bondd[l+1];
          }
          //printf(" deleting bond %i \n",b1);
          nbonds--;
 
          break;
        }
      } //loop k over nbonds

      if (b1>-1) break;
    } //loop j over bonds
  } //loop i over 1+ coordn H

  return;
}

void ICoord::make_bonds()
{
  //printf(" in make_bonds, natoms: %i\n",natoms);
  double MAX_BOND_DIST; 
  nbonds=0;
  for (int i=0;i<natoms;i++)
    for (int j=0;j<i;j++)
    {
       MAX_BOND_DIST = (getR(i) + getR(j))/2;
       if (farBond>1.0) MAX_BOND_DIST *= farBond;
       double d = distance(i,j);
       if (d<MAX_BOND_DIST)
       {
          //printf(" found bond: %i %i dist: %f \n",i+1,j+1,d);
          bonds[nbonds][0]=i;
          bonds[nbonds][1]=j;
          bondd[nbonds]=d;
          nbonds++;
       }
    }

}

void ICoord::coord_num()
{ 
  for (int i=0;i<natoms;i++)
    coordn[i] = 0;
  for (int i=0;i<nbonds;i++)
  {
    coordn[bonds[i][0]]++;
    coordn[bonds[i][1]]++;
  }
}

void ICoord::make_angles()
{
  //include all consecutive connections 
  nangles=0;
  for (int i=0;i<nbonds;i++)
  {
     for (int j=0;j<i;j++)
     {
        if (bonds[i][0]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][1];
          nangles++;
        }
        else if (bonds[i][0]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][0];
          nangles++;
        }
        else if (bonds[i][1]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][1];
          nangles++;
        }
        else if (bonds[i][1]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][0];
          nangles++;
        }
        if (nangles>0)
          anglev[nangles-1]=angle_val(angles[nangles-1][0],angles[nangles-1][1],angles[nangles-1][2]);
     } //loop j
  } //loop i


  return;
}


void ICoord::make_torsions()
{
  int a1,b1,c1,a2,b2,c2;
  bool found;

  ntor = 0;

//  return;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       b1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       b2=angles[j][1];
       c2=angles[j][2];

      // printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,b1,c1,a2,b2,c2);

       if (b1==c2 && b2==c1)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }
       else if (b1==a2 && b2==c1)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       else if (b1==c2 && b2==a1)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }
       else if (b1==a2 && b2==a1)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       if (found && torsions[ntor-1][0] == torsions[ntor-1][2]) { found = false; ntor--; }
       //if (found && bond_exists(torsions[ntor-1][0],torsions[ntor-1][3]))
       if (found && torsions[ntor-1][0] == torsions[ntor-1][3]) { found = false; ntor--; }
       if (found)
       {
         //printf(" made tor: %i %i %i %i \n",torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
         torv[ntor-1]=torsion_val(torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
       }
       //if (found && abs(torv[ntor-1])>180 ) ntor--;
    }
  } 

  return;
}

void ICoord::make_imptor()
{
  int a1,m1,c1,a2,m2,c2;
  bool found;
  nimptor = 0;
  double imptorvt;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       m1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       m2=angles[j][1];
       c2=angles[j][2];

       //printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,m1,c1,a2,m2,c2);

       if (m1==m2)
       {
         if (a1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
         else if (a1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
       } // if m1==m2
       if (found)
       {
//FIX ME
//the following only works when center is 3 coordinate
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1)
             { found = false; nimptor--; }
       }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
//         printf(" imptorv[%i]: %1.4f \n",nimptor,imptorvt);
         if ((abs(imptorvt) < 90.0 || abs(imptorvt) > 120.0) && anumbers[imptor[nimptor-1][2]]==7) { found = false; nimptor--; }
         else if (abs(imptorvt) > 12.0 && abs(imptorvt - 180.) > 12.0) { found = false; nimptor--; }
       }
       if (found) imptorv[nimptor-1] = imptorvt;
    }
  } 

  return;
}

void ICoord::make_imptor_nobonds()
{
  int a1,m1,c1,a2,m2,c2;
  bool found;
  nimptor = 0;
  double imptorvt;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       m1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       m2=angles[j][1];
       c2=angles[j][2];

       //printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,m1,c1,a2,m2,c2);

       if (m1==m2)
       {
         if (a1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
         else if (a1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
       } // if m1==m2
       if (found)
       {
//FIX ME
//the following only works when center is 3 coordinate
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1)
             { found = false; nimptor--; }

         if (anumbers[imptor[nimptor-1][2]]==7)
         {
           //printf(" found imptor centered on nitrogen, deleting \n");
           found = false;
           nimptor--;
         }
       }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
//         printf(" imptorv[%i]: %1.4f \n",nimptor,imptorvt);
//       make all 3 centered atoms planar?
//         printf(" atom: %i has coordn %i \n",imptor[nimptor-1][2],coordn[imptor[nimptor-1][2]]);
         if (coordn[imptor[nimptor-1][2]]!=3)
         {
           found = false;
           nimptor--;
         }
       }
       if (found) imptorv[nimptor-1] = imptorvt;
    }
  } 

  return;
}

double ICoord::torsion_val(int i, int j, int k, int l)
{
  double tval = -999;

  double x1 = coords[3*j+0] - coords[3*i+0];
  double y1 = coords[3*j+1] - coords[3*i+1];
  double z1 = coords[3*j+2] - coords[3*i+2];
  double x2 = coords[3*k+0] - coords[3*j+0];
  double y2 = coords[3*k+1] - coords[3*j+1];
  double z2 = coords[3*k+2] - coords[3*j+2];
  
  double ux1 = y1*z2-z1*y2;
  double uy1 = z1*x2-x1*z2;
  double uz1 = x1*y2-y1*x2;

  double x3 = coords[3*l+0] - coords[3*k+0];
  double y3 = coords[3*l+1] - coords[3*k+1];
  double z3 = coords[3*l+2] - coords[3*k+2];

  double ux2 = z3*y2 - y3*z2;
  double uy2 = x3*z2 - z3*x2;
  double uz2 = y3*x2 - x3*y2;

  double u = (ux1*ux1+uy1*uy1+uz1*uz1)*(ux2*ux2+uy2*uy2+uz2*uz2);

  if (u!=0.0)
  {
     double a = (ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u);
     if (a>1) a=1; else if (a<-1) a=-1;
     tval = acos(a);
     if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
         uz1*(ux2*y2-uy2*x2) < 0.0) tval *=-1;
  }
  else
    tval = 0.;

  if (tval>3.14159) tval-=2*3.14159;
  if (tval<-3.14159) tval+=2*3.14159;

  return tval * 180/3.14;
}

double ICoord::angle_val(int i, int j, int k)
{
   double D1 = distance(i,j);
   double D2 = distance(j,k);
   double D3 = distance(i,k);
   
   double cos = ( D1*D1 + D2*D2 - D3*D3 ) / ( 2*D1*D2);
 
   if (cos > 1) cos = 1;
   if (cos < -1) cos = -1;

  // printf(" cos is: %f \n",cos);
 
   return acos(cos) * 180/3.14159;
}

double ICoord::angle_val_eta2(int i, int j, int k1, int k2)
{
   double D1 = distance(i,j);
//   double D2 = distance(j,k);
//   double D3 = distance(i,k);

  //center of k
   double x1 = (coords[3*k1+0] + coords[3*k2+0])/2.0;
   double x2 = (coords[3*k1+1] + coords[3*k2+1])/2.0;
   double x3 = (coords[3*k1+2] + coords[3*k2+2])/2.0;

   double D2 = sqrt((x1-coords[3*j+0])*(x1-coords[3*j+0])
                   +(x2-coords[3*j+1])*(x2-coords[3*j+1])
                   +(x3-coords[3*j+2])*(x3-coords[3*j+2]));
   double D3 = sqrt((x1-coords[3*i+0])*(x1-coords[3*i+0])
                   +(x2-coords[3*i+1])*(x2-coords[3*i+1])
                   +(x3-coords[3*i+2])*(x3-coords[3*i+2]));

   double cos = ( D1*D1 + D2*D2 - D3*D3 ) / ( 2*D1*D2);
 
   if (cos > 1) cos = 1;
   if (cos < -1) cos = -1;

  // printf(" cos is: %f \n",cos);
 
   return acos(cos) * 180/3.14159;
}

double ICoord::angle_val_eta2_eta2(int i1, int i2, int j, int k1, int k2)
{
//   double D1 = distance(i,j);
//   double D2 = distance(j,k);
//   double D3 = distance(i,k);
 
  //center of i
   double x1 = (coords[3*i1+0] + coords[3*i2+0])/2.0;
   double x2 = (coords[3*i1+1] + coords[3*i2+1])/2.0;
   double x3 = (coords[3*i1+2] + coords[3*i2+2])/2.0;

  //center of k
   double y1 = (coords[3*k1+0] + coords[3*k2+0])/2.0;
   double y2 = (coords[3*k1+1] + coords[3*k2+1])/2.0;
   double y3 = (coords[3*k1+2] + coords[3*k2+2])/2.0;

   double D1 = sqrt((x1-coords[3*j+0])*(x1-coords[3*j+0])
                   +(x2-coords[3*j+1])*(x2-coords[3*j+1])
                   +(x3-coords[3*j+2])*(x3-coords[3*j+2]));
   double D2 = sqrt((y1-coords[3*j+0])*(y1-coords[3*j+0])
                   +(y2-coords[3*j+1])*(y2-coords[3*j+1])
                   +(y3-coords[3*j+2])*(y3-coords[3*j+2]));
   double D3 = sqrt((x1-y1)*(x1-y1)
                   +(x2-y2)*(x2-y2)
                   +(x3-y3)*(x3-y3));

   double cos = ( D1*D1 + D2*D2 - D3*D3 ) / ( 2*D1*D2);
 
   if (cos > 1) cos = 1;
   if (cos < -1) cos = -1;

  // printf(" cos is: %f \n",cos);
 
   return acos(cos) * 180/3.14159;
}

int ICoord::make_nonbond(){

  int n = 0;
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
    {
      bool found = false;
      for (int k=0;k<nbonds;k++)
      {
         if (found) break;
         if ((bonds[k][0]==i && bonds[k][1]==j) ||
             (bonds[k][0]==j && bonds[k][1]==i)) found = true;
      }
      //printf(" checking for pair: %i %i \n",i,j);
      for (int k=0;k<nangles;k++)
      {
        if (found) break;
        //printf(" angle %i bonds: %i %i %i \n",k,angles[k][0],angles[k][1],angles[k][2]);
        if (angles[k][0]==i)
        {
           if (angles[k][1]==j) found = true;
           else if (angles[k][2]==j) found = true;
        }
        else if (angles[k][1]==i)
        {
           if (angles[k][0]==j) found = true;
           else if (angles[k][2]==j) found = true;
        }
        else if (angles[k][2]==i)
        {
           if (angles[k][0]==j) found = true;
           else if (angles[k][1]==j) found = true;
        }
      } // loop k over angles
      if (!found)
      {
        //printf(" not found\n");
        nonbondd[n] = distance(i,j);
        nonbond[n][0] = i;
        nonbond[n][1] = j;
        n++;
      }
    }
  }
  //printf(" n_nonbond: %i \n",n);

  return n;
}

int ICoord::isTM(int a) {

//may later be extended to all 5+ coord types
  int anum;
  if (a>-1)
    anum = anumbers[a];
  else
    return 0;

  int TM = 0;
  if (anum > 1000)
    TM = 2;
  else if (anum > 20)
  {
    if (anum < 31)
      TM = 1;
    else if (38 < anum && anum < 49)
      TM = 1;
    else if (71 < anum && anum < 81)
      TM = 1;
  }
     
  return TM;
}

double ICoord::getR(int i){

  double value = getRa(anumbers,i);

  return value;
}


double ICoord::distance(int i, int j)
{
  //printf("in distance: %i %i\n",i+1,j+1);
  return sqrt((coords[3*i+0]-coords[3*j+0])*(coords[3*i+0]-coords[3*j+0])+
              (coords[3*i+1]-coords[3*j+1])*(coords[3*i+1]-coords[3*j+1])+
              (coords[3*i+2]-coords[3*j+2])*(coords[3*i+2]-coords[3*j+2])); 
}

int ICoord::bond_exists(int b1, int b2) {

   int found = 0;
   if (bond_num(b1,b2)>-1)
     found = 1;
   return found;
}

int ICoord::bond_num(int b1, int b2) {

   int found = -1;

   for (int k1=0;k1<nbonds;k1++)
     if ( (bonds[k1][0] == b1 && bonds[k1][1] == b2)
       || (bonds[k1][1] == b1 && bonds[k1][0] == b2))
     {
       found = k1;
       break;
     }

   return found;
}

int ICoord::hpair(int a1, int a2) {
  if (anumbers[a1]==1 && anumbers[a2]==1)
    return 1;
  else
    return 0;
}

int ICoord::h2count() {

  int count = 0;
  for (int i=0;i<nbonds;i++)
  {
    if (anumbers[bonds[i][0]]==1 && anumbers[bonds[i][1]]==1)
      count++;
  }

  return count;
}




void ICoord::structure_read(string xyzfile)
{  
  //printf("   reading structure \n");  
  
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile)
  {
    printf(" failed to open %s i ICoord \n",xyzfile.c_str());
    exit(-1);
  } 
  
  string line;
  bool success=true;
  success=getline(infile, line);
  if (success)
  {
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
  int natomsa = natoms;
  if (natoms==1) natomsa++;
  printf("   natoms: %2i natomsa: %2i \n",natoms,natomsa);
  if (natoms<1) { printf(" ERROR: natoms must be > 0 \n"); exit(1); }
  
  success=getline(infile, line);
  if (success) comment = line;
  
  anumbers = new int[1+natomsa];
  amasses = new double[1+natomsa];
  anames = new string[1+natomsa];
    
  //cout <<"  -Reading the atomic names...";
  for (int i=0;i<natoms;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i] = tok_line[0];
    anumbers[i] = PTable::atom_number(anames[i]);
    amasses[i] = PTable::atom_mass(anumbers[i]);
  }
  
  infile.close();
  
//  V_profile = new double[1+nnmax];
//  S = new double[1+nnmax];
  
  coords = new double[natomsa*3];
  coords0 = new double[natomsa*3];
   
  //cout <<"  -Reading coordinates...";
 // cout << "Opening the xyz file" << endl;
  infile.open(xyzfile.c_str());
  fflush(stdout);
 // cout << "xyzfile opened" << endl;
  fflush(stdout);
  
  
  success=getline(infile, line);
  success=getline(infile, line);
  for (int j=0;j<natoms;j++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    coords[3*j+0]=atof(tok_line[1].c_str());
    coords[3*j+1]=atof(tok_line[2].c_str());
    coords[3*j+2]=atof(tok_line[3].c_str());
  }
  
  if (natoms!=natomsa)
  {
    coords[3+0] = coords[0]+1.5;
    coords[3+1] = coords[1];
    coords[3+2] = coords[2];
    anames[1] = "X";
    anumbers[1] = 0;
    amasses[1] = 0.;
    natoms = natomsa;
  }
  for (int i=0;i<3*natoms;i++)
     coords0[i] = coords[i];

 // cout << " done" << endl;
  infile.close();
  
 // cout << "Finished reading information from structure file" << endl;
}   

