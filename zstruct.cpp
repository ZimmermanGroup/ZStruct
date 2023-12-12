// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "zstruct.h"
#include "constants.h"

//BUG: 2 coordinate H will break H-transfer (fixed?)
//NOTE: perhaps don't add break move to H in get_h_transfer?

//to do:
// 4. limit max # of reactants (or size)
// 4. (b) set erefs globally (now all vs. elistcmin/single species to initial)
// 5. logistic regression?
// 7. check fragmentation of initial reactants?
// 8. when creating new reactants, label pair_react based on TM in fragment

// 9. "bond" intramolecular fragments




void ZStruct::save_allrxns()
{
  printf("   saving ALLRXNS file. contains ISOMER-TS information \n");

  int* sadd = new int[20];
  int* sbrk = new int[20];


  ofstream arfile;
  string arfile_string = "ALLRXNS";
  arfile.open(arfile_string.c_str());

  int maxadd1 = MAX_ADD+1;
  int maxbrk1 = MAX_BRK+1;

  int nsave = allrxns.size();
  for (int i=0;i<nsave;i++)
  {
    arfile << " " << allrxns[i].id << "  " << allrxns[i].nadd << " " << allrxns[i].nbrks << "   ";

    int natoms1 = allrxns[i].natoms;
    int nadd1 = allrxns[i].nadd;
    int nbrk1 = allrxns[i].nbrks;
//    if (nadd1>maxadd1) nadd1 = maxadd1;
//    if (nbrk1>maxbrk1) nbrk1 = maxbrk1;

    printf(" working on %4i(%4i) with nadd/brk: %i %i \n",i,allrxns[i].id,nadd1,nbrk1);

   //sort adds, brks
    for (int j=0;j<nadd1;j++)
    {
      int a1 = allrxns[i].adds[2*j+0];
      int a2 = allrxns[i].adds[2*j+1];
      if (a1>=natoms1 || a2>=natoms1)
        printf("  ERROR: a1,a2 > natoms: %i %i \n",a1,a2);
      sadd[2*j+0] = a1; sadd[2*j+1] = a2;
      if (allrxns[i].anums[a2]>allrxns[i].anums[a1])
      { sadd[2*j+0] = a2; sadd[2*j+1] = a1; }
    }
    /* printf(" original add:");
    for (int j=0;j<nadd1;j++)
      printf(" %i-%i",sadd[2*j+0],sadd[2*j+1]);
    printf("\n"); */
    for (int j=0;j<nadd1-1;j++)
    {
      if (sadd[2*j+0]<sadd[2*(j+1)+0])
      {
        int t1 = sadd[2*j+0]; int t2 = sadd[2*j+1];
        sadd[2*j+0] = sadd[2*(j+1)+0];
        sadd[2*j+1] = sadd[2*(j+1)+1];
        sadd[2*(j+1)+0] = t1;
        sadd[2*(j+1)+1] = t2;
        j = 0;
      }
    }
    for (int j=0;j<nbrk1;j++)
    {
      int b1 = allrxns[i].brks[2*j+0];
      int b2 = allrxns[i].brks[2*j+1];
      if (b1>=natoms1 || b2>=natoms1)
        printf("  ERROR: b1,b2 > natoms: %i %i \n",b1,b2);
      sbrk[2*j+0] = b1; sbrk[2*j+1] = b2;
      if (allrxns[i].anums[b2]>allrxns[i].anums[b1])
      { sbrk[2*j+0] = b2; sbrk[2*j+1] = b1; }
    }
    /* printf(" original brk:");
    for (int j=0;j<nbrk1;j++)
      printf(" %i-%i",sbrk[2*j+0],sbrk[2*j+1]);
    printf("\n"); */
    for (int j=0;j<nbrk1-1;j++)
    {
      if (sbrk[2*j+0]<sbrk[2*(j+1)+0])
      {
        int t1 = sbrk[2*j+0]; int t2 = sbrk[2*j+1];
        sbrk[2*j+0] = sbrk[2*(j+1)+0];
        sbrk[2*j+1] = sbrk[2*(j+1)+1];
        sbrk[2*(j+1)+0] = t1;
        sbrk[2*(j+1)+1] = t2;
        j = 0;
      }
    }

    /* printf(" rearranged add:");
    for (int j=0;j<nadd1;j++)
      printf(" %i-%i",sadd[2*j+0],sadd[2*j+1]);
    printf("\n");
    printf(" rearranged brk:");
    for (int j=0;j<nbrk1;j++)
      printf(" %i-%i",sbrk[2*j+0],sbrk[2*j+1]);
    printf("\n"); */

    for (int j=0;j<nadd1;j++)
    {
      int a1 = sadd[2*j+0]; int a2 = sadd[2*j+1];
      int an1 = allrxns[i].anums[a1];
      int an2 = allrxns[i].anums[a2];
      int c1 = allrxns[i].coordn[a1];
      int c2 = allrxns[i].coordn[a2];
      double q1 = allrxns[i].q[a1];
      double q2 = allrxns[i].q[a2];
      if (fabs(q1)<0.0000001) q1 = 0.;
      if (fabs(q2)<0.0000001) q2 = 0.;
      printf("  add: %2i %2i  c1,2: %i %i \n",a1,a2,c1,c2);
      arfile << " " << an1 << " " << an2 << " " << c1 << " " << c2 << " " << q1 << " " << q2;
    }
    for (int j=nadd1;j<maxadd1;j++)
      arfile << " 0 0 0 0 0 0";
    arfile << "   ";
    for (int j=0;j<nbrk1;j++)
    {
      int b1 = sbrk[2*j+0]; int b2 = sbrk[2*j+1];
      int an1 = allrxns[i].anums[b1];
      int an2 = allrxns[i].anums[b2];
      int c1 = allrxns[i].coordn[b1];
      int c2 = allrxns[i].coordn[b2];
      double q1 = allrxns[i].q[b1];
      double q2 = allrxns[i].q[b2];
      if (fabs(q1)<0.0000001) q1 = 0.;
      if (fabs(q2)<0.0000001) q2 = 0.;
      printf("  brk: %2i %2i  c1,2: %i %i \n",b1,b2,c1,c2);
      arfile << " " << an1 << " " << an2 << " " << c1 << " " << c2 << " " << q1 << " " << q2;
    }
    for (int j=nbrk1;j<maxbrk1;j++)
      arfile << " 0 0 0 0 0 0";
    arfile << "   ";

    int noverlap = 0;
    for (int j=0;j<nadd1;j++)
    {
      int a1 = sadd[2*j+0]; int a2 = sadd[2*j+1];
      for (int k=0;k<nbrk1;k++)
      {
        int b1 = sbrk[2*k+0]; int b2 = sbrk[2*k+1];
        if (a1==b1 || a1==b2)
          noverlap++;
        if (a2==b1 || a2==b2)
          noverlap++;
      }
    }
    arfile << " " << noverlap;
    arfile << " " << allrxns[i].erxn;
    arfile << "    " << allrxns[i].ets << endl;

  } //loop i over allrxns.size()
  
  arfile.close();

  delete [] sadd;
  delete [] sbrk;

  return;
}

void ZStruct::store_isomers(int id1, int wp1, int nadd1, int nbrks1, int* add1, int* brks1, double* E1, int natoms1, string* anames1, int* anums1, double* xyz1, double* qrc1, int* coordn1)
{
  printf(" adding %4i ",id1);
  irxn1 newrxn;

  newrxn.id = id1;
//need to map pair # to reactant # first
//  newrxn.reactant1 = 
//  newrxn.reactant2 = 

  newrxn.natoms = natoms1;
  newrxn.nadd = nadd1;
  newrxn.nbrks = nbrks1;
  int* add2 = new int[2*nadd1];
  int* brks2 = new int[2*nbrks1];
  for (int i=0;i<2*nadd1;i++) add2[i] = add1[i];
  for (int i=0;i<2*nbrks1;i++) brks2[i] = brks1[i];

  newrxn.adds = add2;
  newrxn.brks = brks2;

  newrxn.estart = E1[0];
  newrxn.ets = E1[2];
  newrxn.erxn = E1[1];

  int* anums2 = new int[natoms1];
  for (int i=0;i<natoms1;i++)
    anums2[i] = anums1[i];
  newrxn.anums = anums2;

  string* anames2 = new string[natoms1];
  for (int i=0;i<natoms1;i++) anames2[i] = anames1[i];
  newrxn.anames = anames2;

  double* q1 = new double[natoms1];
  for (int i=0;i<natoms1;i++) q1[i] = qrc1[i];
  newrxn.q = q1;

  int* coordn2 = new int[natoms1];
  for (int i=0;i<natoms1;i++) coordn2[i] = coordn1[i];
  //for (int i=0;i<natoms1;i++) printf(" cn2: %i",coordn2[i]);
  newrxn.coordn = coordn2;

  for (int i=0;i<2*nadd1;i++)
  if (add2[i]>=natoms1)
    printf(" ERROR: %i > %i \n",add2[i],natoms1);

  int savethis = 1;
  if (natoms1>1000 || nadd1>8 || nbrks1>8)
    savethis = 0;
  for (int i=0;i<2*nadd1;i++)
  if (add1[i]>natoms1)
    savethis = 0;
  for (int i=0;i<2*nbrks1;i++)
  if (brks1[i]>natoms1)
    savethis = 0;
  for (int i=0;i<natoms1;i++)
  if (coordn1[i]>6 || coordn1[i]<0)
    savethis = 0;

  if (savethis)
    allrxns.push_back(newrxn);

  return;
}


void ZStruct::form_reaction_network(int* rxncont, ICoord* icr)
{
  printf("\n in form_reaction_network \n");

  double kT = temperature*(0.6/300.); //need precision
  double Erefts = 50*kT;
  double Erefth = 49*kT;
  double tenkT = 10.*kT;

  printf(" Eref(ts): %3.1f Eref(thermo): %3.1f 10kT: %3.1f \n",Erefts,Erefth,tenkT);

  mark_duplicates(icr);
  write_duplicates(niso);
  get_elistref();

  int* clist = new int[niso];
  for (int i=0;i<niso;i++) clist[i] = -1;

  calc_Ea_Erxn();

 //TS threshold on max Ea
  int nfound = 0;
  for (int i=0;i<niso;i++)
  if (rsaved[i])
  {
    //printf(" iso: %3i \n",i); fflush(stdout);

    double ets = Ea[i];  
    int pass = 0;
    if (ets<Erefts && upair[i]==1)
      pass = 1;

    double erxn = Erxn[i];
    if (pass && erxn>Erefth)
      pass = 0;

    if (pass)
    {
      printf(" ets[%4i]: %12.8f  Initial: %12.8f  Final: %12.8f \n",i,elistts[i],elistref[i],elistp[i]);
      //printf(" ets[%4i]: %4.1f erxn: %5.1f pass: %i \n",i,ets,erxn,pass);
      //if (wp<0) wp = 0;
      //printf(" elistts[i]: %4.1f elistcmin[wp]: %4.1f elistc[i]: %4.1f \n",elistts[i]*627.5,elistcmin[wp]*627.5,elistc[i]*627.5);
    }
    rxncont[i] = pass;

    if (ets<Erefts+5.0 && upair[i]==1 && !pass)
    {
      printf(" ets[%4i]: %12.8f  Initial: %12.8f  Final: %12.8f  near pass \n",i,elistts[i],elistref[i],elistp[i]);
      //printf(" ets[%4i]: %4.1f near pass \n",i,ets);
    }

    Ea[i] = ets;
    if (pass)
      clist[nfound++] = i;
  } //loop i over niso

  int quiet = 1;
 //higher TS/dErxn check
  for (int i=0;i<nfound;i++)
  for (int j=0;j<i;j++)
  {
    int ws1 = clist[i];
    int ws2 = clist[j];
    int wp1 = wpair[ws1];
    int wp2 = wpair[ws2];

    if (wp1==wp2)
    {
      double esh1 = 0.;
      double esh2 = 0.;
      int wsh1 = wshpair[ws1];
      int wsh2 = wshpair[ws2];
      if (wsh1>-1) esh1 = elistr[wshuttle[wshpair[ws1]]];
      if (wsh2>-1) esh2 = elistr[wshuttle[wshpair[ws2]]];

      if (!quiet) printf(" comparing save pair %5i %5i: ",ws1,ws2);

     //previously: kinetics referenced to separated reactants
     //note: kinetics now referenced to lowest complexation species
      double ets1;
      if (wp1>=0 && wsh1==-1)
        ets1 = (elistts[ws1] - elistcmin[wp1])*627.5;
//        ets1 = (elistts[ws1] - elistrp[wp1] - esh1)*627.5;
      else if (wp1>=0 && wsh1>-1)
        ets1 = (elistts[ws1] - elistcminsh[wp1])*627.5;
      else
        ets1 = (elistts[ws1] - elistref[ws1])*627.5; //was elistc
        double ets2;
      if (wp2>=0 && wsh2==-1)
        ets2 = (elistts[ws2] - elistcmin[wp2])*627.5;
//        ets2 = (elistts[ws2] - elistrp[wp2] - esh2)*627.5;
      else if (wp2>=0 && wsh2>-1)
        ets1 = (elistts[ws2] - elistcminsh[wp2])*627.5;
      else
        ets2 = (elistts[ws2] - elistref[ws2])*627.5; //was elistc

      double erxn1;
      if (wp1>=0) 
        erxn1 = (elistp[ws1] - elistrp[wp1] - esh1)*627.5;
      else 
        erxn1 = (elistp[ws1] - elistref[ws1])*627.5; //was elistc
      double erxn2;
      if (wp2>=0)
        erxn2 = (elistp[ws2] - elistrp[wp2] - esh2)*627.5;
      else
        erxn2 = (elistp[ws2] - elistref[ws2])*627.5; //was elistc

      if (erxn1<erxn2-tenkT && ets1<ets2-tenkT)
      {
        if (!quiet) printf(" %5i dominates: %4.1f %4.1f (%5.1f %5.1f) \n",ws1,ets1,ets2,erxn1,erxn2);
//        rxncont[ws2] = 0;
      }
      else if (erxn2<erxn1-tenkT && ets2<ets1-tenkT)
      {
        if (!quiet) printf(" %5i dominates: %4.1f %4.1f (%5.1f %5.1f) \n",ws2,ets1,ets2,erxn1,erxn2);
//        rxncont[ws1] = 0;
      }
      else if (!quiet)
        printf("  pair competes: %4.1f %4.1f (%5.1f %5.1f) \n",ets1,ets2,erxn1,erxn2);
    } //if same reactant pairs
  } //loop i,j over saved


  delete [] clist;

//printf("\n\n ENDING EARLY \n\n");
//exit(1);

  fflush(stdout);

  return;
}

void ZStruct::write_duplicates(int nsave)
{
  printf("   saving DUPLICATES file \n");
  ofstream dfile;
  string dfile_string = "DUPLICATES";
  dfile.open(dfile_string.c_str());

  for (int i=0;i<nsave;i++)
    dfile << upair[i] << endl;

  dfile.close();

  return;
}

int ZStruct::read_duplicates()
{
  string dfilename = "DUPLICATES";
  ifstream output(dfilename.c_str(),ios::in);
  if (!output)
  {
    printf("   didn't find DUPLICATES file \n");
    return 0;
  }

  string line;
  int nfound = 0;
  int isd;
  while (!output.eof())
  {
    getline(output,line);
    //cout << " RR: " << line << endl;
    isd = atoi(line.c_str());
    upair[nfound] = isd;
    nfound++;
    if (nfound>=niso) break;
  }
  output.close();

#if 0
  printf("  printing read-in duplicates \n");
  for (int i=0;i<nfound;i++)
    printf(" %2i \n",upair[i]);
  printf("\n");
#endif

  if (nfound<niso)
  {
    printf("  incomplete duplicates file! found: %4i need: %4i \n",nfound,niso);
#if 0  
    for (int i=0;i<niso;i++) upair[i] = 1;
    for (int i=0;i<niso;i++)
    if (!rsaved[i])
      upair[i] = 0;
#endif
  }

  return nfound;
}

//CPMZ this doesn't check shuttle types
void ZStruct::mark_duplicates(ICoord* icr)
{
  printf("\n in mark_duplicates, niso: %i \n",niso); fflush(stdout);

  if (upair!=NULL) delete [] upair;
  upair = new int[niso+1];
  for (int i=0;i<niso;i++) upair[i] = 1;

  for (int i=0;i<niso;i++)
  if (!rsaved[i])
    upair[i] = 0;

  int dstart = read_duplicates();
  printf("   found %2i rxns in duplicates file \n",dstart);
  if (dstart<0) dstart = 0;
  if (dstart>=niso) return;

  printf(" comparing pairs \n"); fflush(stdout);

 //debug
  //if (niso>30000)
  //  niso = 30000;

  int firstp;
 #pragma omp parallel for private(firstp)
  for (int i=dstart;i<niso;i++)
  {
    //debug
    if (i>50000 && i%10000==0) write_duplicates(i);

    //printf(" working on %4i",i);

    firstp = 1;
    for (int j=0;j<i;j++)
    if (wpair[i]==wpair[j] && upair[i]==1 && upair[j]==1 && wshpair[i] == wshpair[j])
    {
      if (i%100==0 && firstp) { printf(" i: %4i",i); firstp = 0; fflush(stdout); }
      int isD = 0;

      if (natoms[i]>0 && natoms[j]>0)
      if (rsaved[i] && rsaved[j])
      {
        //printf(" (%i-%i)",i,j); fflush(stdout);
        isD = diff_structureiq(natoms[i],anames[i],anumbers[i],xyzp[i],xyzp[j]);
        //isD = diff_structurecq(natoms[i],anames[i],anumbers[i],xyzp[i],xyzp[j]);
        //if (i>9000) { printf(" (%i %i diff: %i) ",i,j,isD); fflush(stdout); }
        if (isD==0)
        {
          int ws = i;
          int wo = j;
          //printf(" duplicates %i %i Ea: %4.1f %4.1f Erxn: %4.1f %4.1f \n",i,j,elistts[i],elistts[j],elistp[i],elistp[j]);
          //print_xyz_gen(natoms[i],anames[i],xyzp[i]);
          //print_xyz_gen(natoms[i],anames[i],xyzp[j]);
          if (elistts[i]<elistts[j])
          {
            ws = j; //mark higher barrier as duplicate
            wo = i;
          }
          upair[ws] = -wo;
        }
      }
      if (natoms[i]==0 || natoms[j]==0)
        upair[i] = 0;
    } //loop j
  } //loop i
  printf("\n\n");

  printf(" comparing to rpairs \n"); fflush(stdout);
//exit(1);

  ICoord icp;
  for (int i=0;i<npair;i++)
  {
    //printf(" checking %i against pair structures \n",i); fflush(stdout);
    int wp = i;
    string nstr = StringTools::int2str(wp,4,"0");
    string pairname = "scratch/rpair"+nstr+".xyz";
    icp.init(pairname);
    int start = pairsnum[2*wp+0];
    int end   = pairsnum[2*wp+1];
    //printf(" start: %i end: %i niso: %i \n",start,end,niso);
    for (int j=start;j<end;j++)
    if (rsaved[j] && natoms[j]>0 && wshpair[j]==-1)
    if (natoms[j]==icp.natoms)
    {
      int isD = diff_structureiq(natoms[j],anames[j],anumbers[j],xyzp[j],icp.coords);
      if (isD==0)
        upair[j] = 100000+i;
    }

   //need to fix memory problem
    //icp.freemem();

  } //loop i over pairs

  printf(" comparing to reactant structures \n"); fflush(stdout);

  for (int i=0;i<nreact;i++)
  {
    printf(" checking against reactant structure %2i \n",i); fflush(stdout);
    int wr = i;
#if 0
    string nstr = StringTools::int2str(wr+1,1,"0");
    string reactname = "react"+nstr+".xyz";
    icp.init(reactname);
#endif
    int start = reactnum[2*wr+0];
    int end   = reactnum[2*wr+1];
    //printf(" start: %i end: %i niso: %i \n",start,end,niso);
    for (int j=start;j<end;j++)
    if (rsaved[j] && natoms[j]>0 && wshpair[j]==-1)
    {
      int isD = diff_structureiq(natoms[j],anames[j],anumbers[j],xyzp[j],icr[wr].coords);
      if (isD==0)
        upair[j] = 100000000+i;
    }

   //need to fix memory problem
    //icp.freemem();
  }


#if 1
  printf(" printing upair:"); fflush(stdout);
  for (int i=0;i<niso;i++)
    printf(" %2i",upair[i]);
  printf("\n");
  fflush(stdout);
#endif
#if 1
  printf(" printing unique:");
  for (int i=0;i<niso;i++)
  if (upair[i]==1)
    printf(" %4i",i);
  printf("\n");
#endif

  return;
}

void ZStruct::create_new_reactants(int* rxncont, ICoord* icr)
{
  printf("\n in create_new reactants, nreact: %i \n",nreact);

  ICoord icr0;
  ICoord icp;

  int nsave = 0;
  for (int i=0;i<niso;i++)
  if (rxncont[i])
    nsave++;

  int nfound = 0;
  int* clist = new int[nsave];
  for (int i=0;i<niso;i++)
  if (rxncont[i])
    clist[nfound++] = i;

  int maxatoms = 0;
  for (int i=0;i<nsave;i++)
  if (natoms[clist[i]]>maxatoms)
    maxatoms = natoms[clist[i]];
  maxatoms *= 2;
  printf("  maxatoms: %i \n",maxatoms);
  double* xyz1 = new double[3*maxatoms];
  int* anumbers1 = new int[maxatoms];
  string* anames1 = new string[maxatoms];
  int* active0 = new int[maxatoms];
  int* active1 = new int[maxatoms];
  for (int i=0;i<maxatoms;i++) active0[i] = 1;

  nfound = 0;
  for (int i=0;i<nsave;i++)
  {
    int ws = clist[i];

    icr0.alloc(natoms[ws]);
    icr0.farBond = 1.1; //CPMZ check this, was using 1.1
    icr0.reset(natoms[ws],anames[ws],anumbers[ws],xyzp[ws]);
    icr0.ic_create();
    icr0.make_frags();
    //icr0.print_xyz();

    int wp0 = wpair[ws];
    int wp = wp0;
    if (wp<0) wp = -wp - 1;
    //printf(" wp: %i ws: %i \n",wp,ws); fflush(stdout);
    //if (pairs==NULL) printf(" pairs is NULL \n"); fflush(stdout);

    int wr1;
    int wr2; 
    if (wp0>-1)
    {
      wr1 = pairs[2*wp+0];
      wr2 = pairs[2*wp+1];
    }
    else
    {
      wr1 = wp;
      wr2 = wp;
    }

    int wsh = wshpair[ws];
    printf("  this product (%2i) from (pair/reactant) %2i reactants: %2i %2i shuttle: %2i \n",ws,wp,wr1,wr2,wsh); fflush(stdout);
    int na1 = natoms[wr1];
    int na2 = natoms[wr2];
    if (wp0<0) na2 = 0;
    int nash = 0;
    if (wsh>-1) nash = natoms[wsh];
    //printf(" na1: %2i na2: %2i nash: %2i total: %2i \n",na1,na2,nash,na1+na2+nash);
#if 1
   //frozen may not carry forward
    for (int j=0;j<na1+na2+nash;j++)
      active0[j] = 1;
    for (int j=0;j<na1;j++)
    if (active[wr1][j]==0)
      active0[j] = 0;
    for (int j=0;j<na2;j++)
    if (active[wr2][j]==0)
      active0[na1+j] = 0;
    for (int j=0;j<nash;j++)
    if (active[wsh][j]==0)
      active0[na1+na2+j] = 0;
#else
    for (int j=0;j<na1;j++)
      active0[j] = active[wr1][j];
    for (int j=0;j<na2;j++)
      active0[na1+j] = active[wr2][j];
    for (int j=0;j<nash;j++)
      active0[na1+na2+j] = active[wsh][j];
#endif
#if 0
    printf(" active0:");
    for (int i=0;i<natoms[ws];i++)
      printf(" %i",active0[i]);
    printf("\n");
#endif

    int nfrags = icr0.nfrags;
    if (nfrags>1)
    {
      int natoms1 = icr0.natoms;
      printf("  splitting up reactant %i into fragments (natoms: %2i) \n",ws,natoms1);
       //icr0.print_xyz();

      for (int j=0;j<nfrags;j++)
      {
        int naf = 0;
        for (int k=0;k<natoms1;k++)
        if (icr0.frags[k]==j)
        {
          xyz1[3*naf+0] = icr0.coords[3*k+0];
          xyz1[3*naf+1] = icr0.coords[3*k+1];
          xyz1[3*naf+2] = icr0.coords[3*k+2];
          anames1[naf]  = icr0.anames[k];
          anumbers1[naf] = icr0.anumbers[k];
          active1[naf]   = active0[k];
          naf++;
        }

        int found = check_against_reactants(naf,anames1,anumbers1,xyz1,icr);
       // printf("  frag %i is duplicate? %i \n",j,found);
        if (!found)
        {
          save_new_reactant(icr,naf,anames1,anumbers1,xyz1,active1);
          nfound++;
        } //if not duplicate
      } //loop j over fragments

    } //if nfrags>1
    else
    {
      printf("  reactant %i is one fragment \n",ws);

      int found = check_against_reactants(natoms[ws],anames[ws],anumbers[ws],xyzp[ws],icr);
      if (!found)
      {
        save_new_reactant(icr,natoms[ws],anames[ws],anumbers[ws],xyzp[ws],active0);
        nfound++;
      }
    }

   //need to fix memory problem
    icr0.freemem();
  } //loop i over nsave

  delete [] clist;
  delete [] anumbers1;
  delete [] anames1;
  delete [] xyz1;
  delete [] active0;
  delete [] active1;

  printf(" done creating new reactants. new: %i total: %i \n",nfound,nreact); fflush(stdout);

  return;
}


void ZStruct::save_new_reactant(ICoord* icr, int natoms1, string* anames1, int* anumbers1, double* xyz1, int* active1)
{
  string nstr = StringTools::int2str(nreact+1,1,"0");
  string reactname = "react"+nstr+".xyz";

  icr[nreact].init(natoms1,anames1,anumbers1,xyz1);
  icr[nreact].print_xyz();
  update_active_mem(nreact+1);
  for (int k=0;k<natoms1;k++)
    active[nreact][k] = active1[k];
#if !DO_NOT_WRITE && !DO_NOT_WRITE_REACTS
  icr[nreact].print_xyz_save(reactname);
  save_frozen(natoms1,active[nreact],nreact+1);
#endif

  nreact++;

  update_react_mem();
  natomsr[nreact-1] = natoms1;

  return;
}

int ZStruct::check_against_reactants(int natoms1, string* anames1, int* anumbers1, double* xyz1, ICoord* icr)
{
  int found = 0;
  for (int wr=0;wr<nreact;wr++)
  {
    //printf(" checking frag against reactant structure %2i \n",wr); fflush(stdout);
    int isD = diff_structureigq(natoms1,icr[wr].natoms,anames1,icr[wr].anames,anumbers1,icr[wr].anumbers,xyz1,icr[wr].coords);
    //int isD = 1;
    if (isD==0) 
    {
      found = 1;
      break;
    }
  } //loop wr over reactants

  return found;
}



int ZStruct::create_pairs_isos(ICoord* icr)
{
  printf("\n\n in create_pairs_isos for %i reactants (niso: %2i) \n",nreact,niso);
  if (nreact<1)
  {
    printf("  not enough isos to create pairs \n");
    return 0;
  }

  int* nisof = new int[nshuttle+2];
  int nisoftot = 0;
  int npfound = 0;

  int npairp = npair;
  int npair0 = nreact*nreact; //was/2+1

  ICoord* icp = new ICoord[npair0];

  update_pairs_mem(npair0,npairp);

#if 1
  printf("\n  before creating pairs \n");
  for (int i=0;i<npair;i++)
    printf("   pair (%i): %i %i range: %i %i \n",i,pairs[2*i+0],pairs[2*i+1],pairsnum[2*i+0],pairsnum[2*i+1]);
  printf("\n");
#endif

  int maxatoms = 0;
  for (int i=0;i<nreact;i++)
  if (icr[i].natoms>maxatoms)
    maxatoms = icr[i].natoms;

  double* xyz1 = new double[6*maxatoms];
  int* anum1 = new int[2*maxatoms];
  string* anames1 = new string[2*maxatoms];
  int* active1 = new int[2*maxatoms];
  double* qrc1 = new double[2*maxatoms];

  for (int i=0;i<nreact;i++)
    printf("  pair_react[%i]: %i \n",i,pair_react[i]);
  for (int i=0;i<nreact;i++)
  for (int j=0;j<=i;j++)
//  if (i!=j)
  if (pair_react[i] || pair_react[j]) //avoid pairing TM's
  {
    int found = 0;
    for (int p1=0;p1<npair;p1++)
    if (pairs[2*p1+0]==i && pairs[2*p1+1]==j)
      found = 1;
    int ool = 0;
    ool = !pair_allowed(i,j,icr);

    if (found)
      printf("   pair %i %i already exists \n",i,j);
    else if (ool)
      printf("   pair %i %i isn't allowed (natoms: %i atomlimit: %i) \n",i,j,icr[i].natoms+icr[j].natoms,atomlimit);
    else
    {
      printf("\n\n pair %i %i is new \n",i,j);
      int nisop = niso;
      int na1 = icr[i].natoms;
      int na2 = icr[j].natoms;
      int q1 = icr[i].q1;
      int q2 = icr[j].q1;
      int natoms1 = na1 + na2;
      printf("  na1: %2i na2: %2i \n",na1,na2);

     // spatially separated fragments
      for (int k=0;k<3*na1;k++)
        xyz1[k]       = icr[i].coords[k];
      for (int k=0;k<3*na2;k++)
        xyz1[3*na1+k] = icr[j].coords[k]+20.;

      for (int k=0;k<na1;k++)
        anames1[k]     = icr[i].anames[k];
      for (int k=0;k<na2;k++)
        anames1[na1+k] = icr[j].anames[k];
      for (int k=0;k<na1;k++)
        anum1[k]     = icr[i].anumbers[k];
      for (int k=0;k<na2;k++)
        anum1[na1+k] = icr[j].anumbers[k];

      for (int k=0;k<na1;k++)
        active1[k]     = active[i][k];
      for (int k=0;k<na2;k++)
        active1[na1+k] = active[j][k];
      for (int k=0;k<na1;k++)
        qrc1[k]     = qr[i][k];
      for (int k=0;k<na2;k++)
        qrc1[na1+k] = qr[j][k];

      printf("  qrc1:");
      for (int k=0;k<natoms1;k++)
        printf(" %4.2f",qrc1[k]);
      printf("\n");

      write_pair_xyz(npair+npfound,natoms1,anames1,xyz1);
      icp[npfound].init(natoms1,anames1,anum1,xyz1);
      icp[npfound].q = qrc1;
      icp[npfound].q1 = q1 + q2;
      icp[npfound].make_frags();
      create_frag_bonds(icp[npfound]);

      printf("\n --- reactions involving pair %2i --- \n",npair+npfound);
      align1.init(na1,icr[i].anames,icr[i].anumbers,icr[i].coords,na2,icr[j].anames,icr[j].anumbers,icr[j].coords);

      for (int k=-1;k<nshuttle;k++)
      {
        nisof[k+1] = create_isos(icp[npfound],2,niso,active1,k); //niso incremented within
        nisoftot += nisof[k+1];
        printf(" after create_isos, wsh: %i nisof: %i \n",k,nisof[k+1]);
      }

      align1.freemem();
      printf("\n --- done with (%i) reactions involving pair %2i --- \n\n",nisoftot,npair+npfound);
 
      update_wpair_mem(nisop,nisoftot);
      printf(" nisop+nisoftot: %3i niso: %3i \n",nisop+nisoftot,niso);
      for (int k=nisop;k<niso;k++)
        wpair[k] = npair+npfound;
      int shstart = nisop;
      for (int k=-1;k<nshuttle;k++)
      {
        int shend = shstart+nisof[k+1];
        printf("   wsh: %i start: %i end: %i \n",k,shstart,shend);
        for (int l=shstart;l<shend;l++)
          wshpair[l] = k;
        shstart = shend;
      }
      pairs[2*(npair+npfound)+0] = i;
      pairs[2*(npair+npfound)+1] = j;
      pairsnum[2*(npair+npfound)+0] = nisop;
      pairsnum[2*(npair+npfound)+1] = niso;
      npfound++;
    }//if !found
  } //if allowed

#if 1
  printf("\n  after creating pairs \n");
  for (int i=0;i<npair+npfound;i++)
    printf("   pair (%i): %i %i range: %i %i \n",i,pairs[2*i+0],pairs[2*i+1],pairsnum[2*i+0],pairsnum[2*i+1]);
  printf("\n");
#endif
#if 0
  printf(" wpairs:");
  for (int i=0;i<niso;i++)
    printf(" %2i",wpair[i]);
  printf("\n");
  printf(" wshpairs:");
  for (int i=0;i<niso;i++)
    printf(" %2i",wshpair[i]);
  printf("\n");
#endif

  delete [] nisof;
  delete [] xyz1;
  delete [] anames1;
  delete [] anum1;
  delete [] active1;
  delete [] qrc1;
  for (int i=0;i<npfound;i++)
    icp[i].freemem();
  delete [] icp;

//printf("\n\n ENDING EARLY \n");
//exit(1);

  npair = npfound + npairp;
  return nisoftot;
}


int ZStruct::create_react_isos(ICoord* icr)
{
  int* nisof = new int[nshuttle+1];
  int nisoftot = 0;

#if 1
  printf("\n before create_react_isos \n");
  for (int i=0;i<nreact;i++)
    printf("  reactnum for %i: %i %i \n",i,reactnum[2*i+0],reactnum[2*i+1]);
#endif

  int natomsmax = 0;
  for (int i=0;i<nreact;i++)
  if (natomsr[i]>natomsmax)
    natomsmax = natomsr[i];
  natomsmax *= 2;
  double* qrc1 = new double[natomsmax];

  nfragb = 0;
  for (int i=0;i<nreact;i++)
  if (reactrun[i]==0)
  {
    printf("\n\n");
    printf("\n --- doing reactant %i (niso: %i) --- \n",i+1,niso);
    int nisop = niso;
    reactnum[2*i+0] = niso;

    for (int j=0;j<natomsr[i];j++)
      qrc1[j] = qr[i][j];
    icr[i].q = qrc1;

    printf("  qrc1:");
    for (int j=0;j<natomsr[i];j++)
      printf(" %4.2f",qrc1[j]);
    printf("\n");

    for (int j=-1;j<nshuttle;j++)
    {
      nisof[j+1] = create_isos(icr[i],1,niso,active[i],j);
      nisoftot += nisof[j+1];
    }
    reactnum[2*i+1] = niso;

    update_wpair_mem(nisop,niso);
    int shstart = nisop;
    for (int k=-1;k<nshuttle;k++)
    {
      int shend = shstart+nisof[k+1];
      printf("   wsh: %i start: %i end: %i \n",k,shstart,shend);
      for (int l=shstart;l<shend;l++)
        wshpair[l] = k;
      shstart = shend;
    }
    for (int k=nisop;k<niso;k++)
      wpair[k] = -i-1;
    nisop = niso;
  }
  else
    printf(" reactant %i already done \n",i+1);

  printf("\n after create_react_isos \n");
  for (int i=0;i<nreact;i++)
    printf("  reactnum for %i: %i %i \n",i+1,reactnum[2*i+0],reactnum[2*i+1]);

  delete [] qrc1;
  delete [] nisof;

  return nisoftot;
}

void ZStruct::do_one_step(ICoord* icr)
{
  update_react_mem();

  int nisof = 0;
  if (nreact>0)
  {
    nisof = create_react_isos(icr);
#if DO_PAIR_ISOS
    nisof += create_pairs_isos(icr);
#endif
  }
  else
  {
    printf(" ERROR: no reactants! \n");
    exit(1);
  }
  printf("\n found %3i total isomers, %3i new \n",niso,nisof);

  dft_para(nreact,icr);
  fill_elistrp();

#if USE_NBO
  get_nbo_reactants(icr);
#endif

  save_gsm(0);
  save_shuttles(0);
  gsm_para(niso-nisof,niso);
  save_gsm(1);
  save_shuttles(1);

  return;
}

void ZStruct::get_one_nbo(int a1, int type, NBO* ntest1, Mopac& mopt1)
{
  if (type!=1 && type!=2)
  {
    printf("  ERROR: type should be 1 (reactant) or 2 (product) \n");
    exit(1);
  }

#if !USE_MOPAC
  int natomsmax = 0;
  for (int i=0;i<niso;i++)
  if (natomsmax<natoms[i])
    natomsmax = natoms[i];

  DFT dft1;
  dft1.alloc(natomsmax);
#endif

  string filename_in, filename_out;
#if USE_MOPAC
  string nstr = StringTools::int2str(a1,4,"0");
  if (type==1)
  {
    filename_in = "scratch/startmop"+nstr;
    filename_out = "scratch/startmop"+nstr+".out";
  }
  else if (type==2)
  {
    filename_in = "scratch/prodmop"+nstr;
    filename_out = "scratch/prodmop"+nstr+".out";
  }
#elif 1
  string nstr = StringTools::int2str(a1,4,"0");
  if (type==1)
  {
    filename_in = "scratch/startdft"+nstr;
    filename_out = "scratch/startdft"+nstr+".out";
  }
  else if (type==2)
  {
    filename_in = "scratch/proddft"+nstr;
    filename_out = "scratch/proddft"+nstr+".out";
  }
#else
  string nstr = StringTools::int2str(a1,4,"0");
  if (type==1)
  {
    filename_in = "scratch/*/startdft"+nstr;
    filename_out = "scratch/*/startdft"+nstr+".out";
  }
  else if (type==2)
  {
    filename_in = "scratch/*/proddft"+nstr;
    filename_out = "scratch/*/proddft"+nstr+".out";
  }
#endif

  if (type==1)
    ntest1[a1].init(natoms[a1],anumbers[a1],anames[a1],xyzr[a1]);
  else if (type==2)
    ntest1[a1].init(natoms[a1],anumbers[a1],anames[a1],xyzp[a1]);

  struct stat sts;
  if (stat(filename_out.c_str(), &sts) != -1)
  {
    //printf("   already calculated product %2i NBO \n",i+1);
  }
  else
  {
#if USE_MOPAC
    //printf("   need to get NBO for product %2i \n",a1+1);
    if (type==1)
      mopt1.reset(natoms[a1],anumbers[a1],anames[a1],xyzr[a1]);
    else if (type==2)
      mopt1.reset(natoms[a1],anumbers[a1],anames[a1],xyzp[a1]);
    double energy = mopt1.sp_energy(filename_in);
#else
    printf(" DFT preparing structure for sp ");
    if (type==1)
      dft1.reset(natoms[a1],anumbers[a1],anames[a1],xyzr[a1]);
    else if (type==2)
      dft1.reset(natoms[a1],anumbers[a1],anames[a1],xyzp[a1]);
#if !DO_NOT_WRITE
    dft1.sp_dnr(filename_in);
#endif
#endif
  }

  int nbotype = USE_NBO;
  if (nbotype==2)
    return;

  int pass = ntest1[a1].read_nbo_file(filename_out);
  //ntest1[a1].print_nbo();
  if (!pass) printf("  reading NBO file %s failed \n",filename_out.c_str()); 

#if !USE_MOPAC
  dft1.freemem();
#endif

  return;
}

void ZStruct::get_nbo_reactions(ICoord* icr)
{
  printf("\n\n get_nbo_reactions:  NOT YET IMPLEMENTED \n");
  printf("  niso: %2i \n",niso);

  int nbotype = USE_NBO;

  int natomsmax = 10;
  for (int i=0;i<niso;i++)
  if (natoms[i]>natomsmax)
    natomsmax = natoms[i];

  Mopac mopt;
  mopt.alloc(natomsmax);

 //reactant complexes
  NBO* ntest0 = new NBO[niso];
  for (int i=0;i<niso;i++)
  //for (int i=508;i<509;i++)
  if (rsaved[i])
  {
    //printf(" working on reactant %2i \n",i+1);
    get_one_nbo(i,1,ntest0,mopt);
  }

 //products
  NBO* ntest = new NBO[niso];
  for (int i=0;i<niso;i++)
  //for (int i=508;i<509;i++)
  if (rsaved[i])
  {
    //printf(" working on product %2i \n",i+1);
    get_one_nbo(i,2,ntest,mopt);
  }


//  if (nbotype==2) return;

  //printf("\n\n ENDING EARLY \n");
  //exit(1);  

  int* oadd = new int[niso];
  int* obrk = new int[niso];
  double* addq = new double[8*niso];
  double* brkq = new double[8*niso];
  double* addspr = new double[8*niso];
  double* brkspr = new double[8*niso];
  double* addpol = new double[4*niso];
  double* addpol2 = new double[4*niso];
  double* brkpol = new double[4*niso];
  double* addeig = new double[4*niso];
  double* addeig2 = new double[4*niso];
  double* brkeig = new double[4*niso];
  for (int i=0;i<niso;i++) oadd[i] = obrk[i] = 0;
  for (int i=0;i<4*niso;i++) addpol[i] = brkpol[i] = 0.;
  for (int i=0;i<4*niso;i++) addeig[i] = brkeig[i] = 0.;
  for (int i=0;i<8*niso;i++) addq[i] = brkq[i] = 0.;
  for (int i=0;i<8*niso;i++) addspr[i] = brkspr[i] = 0.;


#if 0
  int i1 = 318;
  //for (int i=i1;i<i1+1;i++)
  for (int i=0;i<niso;i++)
  {
    printf("\n\n Reaction #%4i Ea: %5.1f \n",i,Ea[i]);
   //print reactions
    ntest0[i].print_nbo();
    ntest[i].print_nbo();
    printf("  broken: \n");
    ntest0[i].compare_nbo(ntest[i],"none",0);
    printf("  formed: \n");
    ntest[i].compare_nbo(ntest0[i],"none",0);
  }
#endif

  if (nbotype<2)
  for (int i=0;i<niso;i++)
  if (rsaved[i])
  {
#if 0
    printf("\n\n Reaction #%4i Ea: %5.1f \n",i,Ea[i]); fflush(stdout);
   //print reactions
    ntest0[i].print_nbo();
    ntest[i].print_nbo();
    printf("  broken: \n");
    ntest0[i].compare_nbo(ntest[i],"none",0);
    printf("  formed: \n");
    ntest[i].compare_nbo(ntest0[i],"none",0);
#else
   //quiet mode
    ntest0[i].compare_nbo(ntest[i],"none",1);
    ntest[i].compare_nbo(ntest0[i],"none",1);
#if !USE_MOPAC
    ntest0[i].blist = new int[10]; for (int j=0;j<10;j++) ntest0[i].blist[j] = 0;
    ntest[i].blist = new int[10];  for (int j=0;j<10;j++) ntest[i].blist[j] = 0;
#endif
#endif

#if 0
    for (int j=0;j<min(ntest0[i].nb,4);j++)
    {
      int b1 = ntest0[i].blist[j];
      int bt1 = ntest0[i].bmo_atoms[2*b1+0]-1;
      int bt2 = ntest0[i].bmo_atoms[2*b1+1]-1;
      if (bt2<0) bt2 = bt1;
      obrk[i]++;

      brkq[8*i+2*j+0] = ntest0[i].q[bt1];
      brkq[8*i+2*j+1] = ntest0[i].q[bt2];
      brkspr[8*i+2*j+0] = ntest0[i].spratio[bt1];
      brkspr[8*i+2*j+1] = ntest0[i].spratio[bt2];
      sort_2(&brkq[8*i+2*j]);
      sort_2(&brkspr[8*i+2*j]);

      brkpol[4*i+j] = ntest0[i].bmo_polar[b1];
      //brkpol2[4*i+j] = brkq[8*i+2*j+0] - brkq[8*i+2*j+1];
      brkeig[4*i+j] = ntest0[i].mo_occ[b1]/10.; //set to ~1.
      //printf(" orbital #: %2i brkpol: %6.3f brkpol2: %6.3f brkeig: %5.3f \n",b1,brkpol[4*i+j],brkpol[4*i+j],brkeig[4*i+j]);
    }
    for (int j=0;j<min(ntest[i].nb,4);j++)
    {
      int a1 = ntest[i].blist[j];
      int at1 = ntest[i].bmo_atoms[2*a1+0]-1;
      int at2 = ntest[i].bmo_atoms[2*a1+1]-1;
      if (at2<0) at2 = at1;
     //get orbitals on atoms at1 and at2?
      oadd[i]++;

      addq[8*i+2*j+0] = ntest0[i].q[at1];
      addq[8*i+2*j+1] =	ntest0[i].q[at2];
      addspr[8*i+2*j+0] = ntest0[i].spratio[at1];
      addspr[8*i+2*j+1] = ntest0[i].spratio[at2];
      sort_2(&addq[8*i+2*j]);
      sort_2(&addspr[8*i+2*j]);

      addpol[4*i+j] = ntest[i].bmo_polar[a1];
      addpol2[4*i+j] = addq[8*i+2*j+0] - addq[8*i+2*j+1];
      addeig[4*i+j] = ntest[i].mo_occ[a1]/10.; //set to ~1.
      //printf(" orbital #: %2i addpol: %6.3f addpol2: %6.3f addeig: %5.3f \n",a1,addpol[4*i+j],addpol2[4*i+j],addeig[4*i+j]);
    }
#endif

   //putting into standard order
    if (ntest0[i].nb==2) sort_2p(&brkpol[4*i],&brkeig[4*i]);
    if (ntest[i].nb==2)  sort_2p(&addpol2[4*i],&addeig[4*i]);
    if (ntest0[i].nb==3) sort_3p(&brkpol[4*i],&brkeig[4*i]);
    if (ntest[i].nb==3)  sort_3p(&addpol2[4*i],&addeig[4*i]);
    if (ntest0[i].nb==4) sort_4p(&brkpol[4*i],&brkeig[4*i]);
    if (ntest[i].nb==4)  sort_4p(&addpol2[4*i],&addeig[4*i]);
  }

//  exit(1);

  double elim = 100;
  printf("\n\n comparing reactions with barriers less than %4.1f \n",elim);
  fflush(stdout);

#if 0
  for (int i=0;i<3;i++)
  for (int j=0;j<i;j++)
  if (Ea[i] < elim && Ea[j] < elim)
  {
    printf("  %4i vs. %4i \n",i,j);
    printf("   oadd: %i/%i obrk: %i/%i addpol:",oadd[i],oadd[j],obrk[i],obrk[j]);
    for (int k=0;k<oadd[i];k++)
      printf(" %5.3f",addpol[4*i+k]);
    printf("/");
    for (int k=0;k<oadd[j];k++)
      printf(" %5.3f",addpol[4*j+k]);
    printf(" brkpol:");
    for (int k=0;k<obrk[i];k++)
      printf(" %5.3f",brkpol[4*i+k]);
    printf("/");
    for (int k=0;k<obrk[j];k++)
      printf(" %5.3f",brkpol[4*j+k]);
    printf(" addeig:");
    for (int k=0;k<oadd[i];k++)
      printf(" %5.3f",addeig[4*i+k]);
    printf("/");
    for (int k=0;k<oadd[j];k++)
      printf(" %5.3f",addeig[4*j+k]);
    printf(" brkeig:");
    for (int k=0;k<obrk[i];k++)
      printf(" %5.3f",brkeig[4*i+k]);
    printf("/");
    for (int k=0;k<obrk[j];k++)
      printf(" %5.3f",brkeig[4*j+k]);
    printf(" Ea: %5.2f/%5.2f \n",Ea[i],Ea[j]);
//    printf("\n");
  }
#endif

  int* ids = new int[niso];
  for (int i=0;i<niso;i++) ids[i] = i;
  for (int i=0;i<niso;i++)
  if (rsaved[i]<1)
    ids[i] = niso+1000000;
//  double* dist = new double[niso*niso];
#if 0
  double alpha = 5.;
  for (int i=0;i<niso;i++)
  for (int j=0;j<niso;j++)
  {
    double d0 = 0.;
    if (oadd[i]==oadd[j] && obrk[i]==obrk[j])
    {
      for (int k=0;k<min(oadd[i],4);k++)
      {
        d0 += fabs(addpol2[4*i+k] - addpol2[4*j+k]);
        //d0 += fabs(addeig2[4*i+k] - addeig2[4*j+k]);
        d0 += fabs(addspr[8*i+2*k+0] - addspr[8*j+2*k+0]);
        d0 += fabs(addspr[8*i+2*k+1] - addspr[8*j+2*k+1]);
      }
      for (int k=0;k<min(obrk[i],4);k++)
      {
        d0 += fabs(brkpol[4*i+k] - brkpol[4*j+k]);
        d0 += fabs(brkeig[4*i+k] - brkeig[4*j+k]);
        d0 += fabs(brkspr[8*i+2*k+0] - brkspr[8*j+2*k+0]);
        d0 += fabs(brkspr[8*i+2*k+1] - brkspr[8*j+2*k+1]);
      }
    }
    else d0 = 50.;
    dist[i*niso+j] = d0;
    //double d1 = exp(-d0*d0/alpha);
   // printf(" %4i/%4i d: %8.5f %8.5f \n",i,j,d0,d1);
  }
#endif

 //binary unique reactions data
  int* urxn1 = new int[niso];
  for (int i=0;i<niso;i++) urxn1[i] = 0;
  for (int i=0;i<niso;i++)
  if (upair[i]==1)
    urxn1[i] = 1;
  //for (int i=0;i<niso;i++)
  //if (upair[i]!=1)
  //  ids[i] += 1000000;

  rxns2.add_nbo_data(niso,ntest0,ntest,ids);
  rxns2.set_unique(niso,urxn1,ids);
  rxns2.temperature = temperature;


 //allocates space for regressions
  rxns2.init_ga(1);
  rxns2.init_ga_2();
  double err = rxns2.create_regression_nbo(2);
  rxns2.save_xy_data("xydata");

 exit(1);

 //SIGMOID setup
  double Ea_sig = 20.;
  double alpha_sig = 10.;
  double* expEa = new double[niso];
  for (int i=0;i<niso;i++) expEa[i] = sigmoid(Ea[i],Ea_sig,alpha_sig);
  for (int i=0;i<niso;i++) if (Ea[i]>111.1) expEa[i] = 111.1;
//  save_gephi(niso,ids,dist,expEa,5.0,0.25);



#if 0
  printf(" Printing Kernel matrix! \n");
  for (int i=0;i<niso;i++)
  {
    for (int j=0;j<niso;j++)
      printf(" %12.8f",dist[i*niso+j]);
    printf("\n");
  }
#endif




  int nftr = 28;
  double* X = new double[niso*nftr];
  for (int i=0;i<niso*nftr;i++) X[i] = 0.;
  for (int i=0;i<niso;i++)
  {
    for (int j=0;j<4;j++)
      X[i*nftr+j] = addpol[4*i+j];
    for (int j=0;j<4;j++)
      X[i*nftr+4+j] = brkpol[4*i+j];
    for (int j=0;j<8;j++)
      X[i*nftr+8+j] = addspr[8*i+j];
    for (int j=0;j<8;j++)
      X[i*nftr+16+j] = brkspr[8*i+j];
    for (int j=0;j<4;j++)
      X[i*nftr+24+j] = brkpol[4*i+j];
  }
  KNNR knnr1;
  knnr1.init();
  //knnr1.reset_active();
  knnr1.ids = ids;

  for (int i=0;i<niso;i++) if (Ea[i]>111.1) expEa[i] = 0.;
  knnr1.load_values(niso,nftr,X,expEa);
  double KNNR_PTHRESH = 0.5;
  knnr1.pthresh = KNNR_PTHRESH;

  int KNNR_K = 4;
  double error = knnr1.test_points(KNNR_K);
  printf(" error: %8.5f \n",error);

  printf("\n printing false negatives \n");
  for (int i=0;i<niso;i++)
  {
    double val = knnr1.errlist[i] + expEa[i];
    if (val<KNNR_PTHRESH && expEa[i]>KNNR_PTHRESH) 
    {
      printf("  pt: %4i id: %4i val: %3.2f actual: %3.2f \n",i,ids[i],val,expEa[i]);
      for (int l=0;l<KNNR_K;l++)
      {
        int index = knnr1.knnlist[i*KNNR_K+l];
        for (int j=0;j<nftr;j++)
          printf(" %4.1f",X[index*nftr+j]*10.);
        printf("  %4.3f   id: %3i/%3i",expEa[index],index,ids[index]);
        printf("  dist: %4.3f \n",knnr1.get_distance(&X[i*nftr],&X[index*nftr]));
      }
    }
  }



  delete [] oadd;
  delete [] obrk;
  delete [] addq;
  delete [] brkq;
  delete [] addspr;
  delete [] brkspr;
  delete [] addpol;
  delete [] addpol2;
  delete [] brkpol;
  delete [] addeig;
  delete [] addeig2;
  delete [] brkeig;

//  delete [] dist;
  delete [] ntest;

  printf(" done with get_nbo_reactions \n");
  return;
}

void ZStruct::get_nbo_reactants(ICoord* icr)
{
//  return;

  printf("\n\n\n in get_nbo for %2i reactants \n",nreact);

  string filename;

  if (nbo_reacts!=NULL) delete [] nbo_reacts;
  nbo_reacts = new NBO[nreact];

  for (int i=0;i<nreact;i++)
  {
    printf(" working on reactant %2i \n",i+1);
#if USE_MOPAC
    string nstr = StringTools::int2str(i,4,"0");
    filename = "scratch/reactmop"+nstr+".out";
#else
    string nstr = StringTools::int2str(i,4,"0");
    filename = "scratch/reactdft"+nstr+".out";
#endif
    nbo_reacts[i].init(icr[i].natoms,icr[i].anumbers,icr[i].anames,icr[i].coords);
    nbo_reacts[i].read_nbo_file(filename);
    nbo_reacts[i].print_nbo();
  }

#if ISOMER_DB && USE_MOPAC
  if (qr!=NULL) delete [] qr;
  qr = new double*[nreact];

  for (int i=0;i<nreact;i++)
  {
    printf(" putting charge into q \n");
    int natoms1 = icr[i].natoms;
    qr[i] = new double[natoms1];
    for (int j=0;j<natoms1;j++)
      qr[i][j] = nbo_reacts[i].q[j];
  }
#endif

#if 0
  printf("  comparing react1 to react2 \n");
  filename = "m1vs2.out";
  ntest[0].compare_nbo(ntest[1],filename);
  printf("  comparing react2 to react1 \n");
  filename = "m2vs1.out";
  ntest[1].compare_nbo(ntest[0],filename);
#endif

  return;
}


void ZStruct::go_zstruct(int nsteps)
{
  printf("  **** test function go_zstruct **** \n\n");
  printf(" nsteps: %i \n",nsteps);


  init();
  int nelemf = get_limits();
  temperature = read_temperature();
  if (temperature < 0.) temperature = 298.15;
  printf(" T set to: %4.2f \n\n",temperature);

  printf("\n MAX_ADD: %i MAX_BRK: %i MAX_CHG: %i MIN_CHG: %i \n",MAX_ADD,MAX_BRK,MAX_CHG,MIN_CHG);
#if DO_NOT_WRITE
  printf(" DO_NOT_WRITE is ON \n");
#endif
#if USE_H_TRANSFER
  printf(" USE_H_TRANSFER is ON \n");
#else
  printf(" USE_H_TRANSFER is OFF \n");
#endif
  printf(" RINGSHUTTLE is %i \n",RINGSHUTTLE);
#if DO_PAIR_ISOS
  printf(" DO_PAIR_ISOS is ON \n");
#else
  printf(" DO_PAIR_ISOS is OFF \n");
#endif
#if TM_PAIRS
  printf(" TM_PAIRS is ON \n");
#else
  printf(" TM_PAIRS is OFF \n");
#endif
  printf("\n");

  int maxr = 500;
  ICoord* icr = new ICoord[maxr];
  nreact = read_reactants(icr);

  nshuttle = read_shuttles();
  printf(" found %i shuttles \n",nshuttle);
  icshuttle = new ICoord[nshuttle+1];
  for (int i=0;i<nshuttle;i++)
  {
    int wr = wshuttle[i];
    icshuttle[i].init(icr[wr].natoms,icr[wr].anames,icr[wr].anumbers,icr[wr].coords);
  }


 //also allocates "active" 
  int nfrz = get_frozen(nreact);
  //handle "X"
  for (int i=0;i<nreact;i++)
  for (int j=0;j<icr[i].natoms;j++)
  if (icr[i].anames[j]=="X")
    active[i][j] = 0;

  dft_para(nreact,icr); //update_elistrp done in read_gsm_data

#if ISOMER_DB
  get_nbo_reactants(icr);
#endif

  int dataf = read_gsm_data();
  printf(" found %2i potential GSM files \n\n",dataf); fflush(stdout);

#if ISOMER_DB
  save_allrxns();
//  exit(1);
#endif

  //rxns2.print_reactions();
#if USE_GA
  if (dataf>0) form_rxn_ga();
#endif
  niso += dataf;

#if USE_DT
  rxns2.temperature = temperature;
  rxns2.make_decision_tree();
#endif

#if USE_GA
  rxns2.element_analysis(-1);
#endif
#if TEST_GA && USE_GA

printf("\n ENDING EARLY \n");
exit(1);
#if 0
  printf("\n\n ******** TEST MODE RESETTING PAIRS ******** \n\n");
  niso = 155;
  for (int i=3;i<npair;i++)
    pairs[2*i+0] = pairs[2*i+1] = -1;
  for (int i=2;i<nreact;i++)
    reactrun[i] = 0;
  npair = 3;
#elif 0
  printf("\n\n ******** TEST MODE RESETTING PAIRS ******** \n\n");
  niso = 0;
  for (int i=0;i<npair;i++)
    pairs[2*i+0] = pairs[2*i+1] = -1;
  for (int i=0;i<nreact;i++)
    reactrun[i] = 0;
  //npair = 0;
#endif
#endif


  int nbotype = USE_NBO;


  int* rxncont = NULL;
  for (int n=0;n<nsteps;n++)
  {
    printf("\n\n  Starting step #%i \n",n+1);
    if (nbotype<2)
      do_one_step(icr);

    dataf = read_gsm_data();
    printf(" after read, niso: %i dataf: %i \n",niso,dataf);

#if USE_GA
    if (dataf>0) form_rxn_ga();
//printf(" ENDING EARLY \n");
//exit(1);
#endif
    niso = dataf;

    if (rxncont!=NULL) delete [] rxncont;
    rxncont = new int[niso];
    for (int i=0;i<niso;i++) rxncont[i] = 0;

    form_reaction_network(rxncont,icr);
#if 0
    create_new_reactants(rxncont,icr);
#endif

  } //loop n over nsteps

#if 0
  dataf = read_gsm_data();
  printf(" found %2i potential GSM files \n",dataf);

  //
#endif

  if (nbotype>0)
  {
    get_nbo_reactions(icr);

    printf("\n\n ENDING EARLY \n");
    exit(1);
  }



  printf("\n  **** done with test go_z **** \n");


  delete [] rxncont;
#if 0
  for (int i=0;i<nreact;i++)
    icr[i].freemem();
  delete [] icr;
#endif

  return;
}

int ZStruct::ml_eval(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1)
{
  double prob0 = rxns2.tree_screen_1(ic1.natoms,ic1.anames,ic1.anumbers,ic1.coords,nadd1,add1,nbrks1,brks1,ic1.q);
#if USE_DT
  printf(" prob0: %4.3f \n",prob0);
#endif

#if USE_GA
  return ga_eval(ic1,nadd1,add1,nbrks1,brks1);
#endif

  int pass = 0;
  if (prob0>0.5)
    pass = 1;
 
  return pass;
}

int ZStruct::ga_eval(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1)
{
  if (ga_avail==0) return 1;

  int pass = 1;
  //printf("  in ga_eval nadd: %2i nbrks: %2i",nadd1,nbrks1); fflush(stdout);

  double prob = rxns2.eval_prob(ic1.natoms,ic1.anames,ic1.anumbers,ic1.coords,nadd1,add1,nbrks1,brks1,ic1.q);

  //if (prob==-999.) printf(" no prediction available \n");
  if (prob < pthresh) pass = 0;
  if (prob!=-999.)
  {
#if TEST_GA
    double ets;
    double esh = 0.; 
    int i = niso;
    int wp = wpair[niso];
    int wsh = wshpair[i];
    if (wsh>-1) esh = elistr[wshuttle[wshpair[i]]];
    //printf(" esh: %8.6f ",esh);
    //printf(" wp: %i ",wp);
//    if (wp>=0)
//      ets = (elistts[i] - elistrp[wp] - esh)*627.5;
    if (wp>=0 && wsh==-1)
      ets = (elistts[i] - elistcmin[wp])*627.5;
    else if (wp>=0)
      ets = (elistts[i] - elistcminsh[wp])*627.5;
    else
      ets = (elistts[i] - elistc[i])*627.5;
    double simfr = -1.;
    int wsdb = rxns2.find_ts(i);
    if (wsdb>-1) simfr = rxns2.similar[wsdb];
    printf("    iso: %3i prob: %6.3f pass? %i Ea: %3.1f Sim: %4.3f \n",niso,prob,pass,ets,simfr);
//    printf("    iso: %3i prob: %6.3f pass? %i Ea: %3.1f \n",niso,prob,pass,ets);
#else
    printf("    iso: %3i prob: %6.3f pass? %i \n",niso,prob,pass);
#endif
  }
  else
  {
    pass = 1;
    printf("    iso: %3i prob:   N/A  pass? %i \n",niso,pass);
  }


 //CPMZ debug
#if ALL_PASS
  pass = 1;
#endif

  //printf("  done ga_eval \n"); fflush(stdout);

  return pass;
}

void ZStruct::form_rxn_ga()
{
  rxns1.temperature = temperature;
  rxns2.temperature = temperature;

#if 0
  printf(" forming regression for rxns1 \n");
  rxns1.create_regression(0);
#elif 0
  printf("\n forming regression for rxns2 \n");
  rxns2.create_regression(0);
#endif

  int ga_pop_1 = GA_POP;
  int ga_steps_1 = GA_STEPS;
  int nga = 0;
  if (!ga_avail) nga = read_ga();
  if (nga==ga_pop_1)
  {
    printf(" found GA! \n");
    ga_avail = 1;
//    ga_steps_1 = 1; //CPMZ note
    rxns2.init_ga(nga);
    rxns2.init_ga_2();
  }
  else
  {
    //printf(" couldn't find ga save file \n\n");
    ga_avail = 0;
  }

  printf("\n\n Starting GA  \n");
  if (!ga_avail)
    rxns2.init_ga_features(ga_pop_1);
  rxns2.ga_run(ga_steps_1);

 //not implemented
// rxns2.element_analysis(); 


#if 0
 //debug
  int wi = 3; 
  printf("\n\n evaluating probability for %i:",wi);
  int success = read_ISOMER(wi,nadd,add,nbrks,brks);
  int nnodes = read_string(wi,energies,xyz);
  double prob = rxns2.eval_prob(ic1.natoms,ic1.anames,ic1.anumbers,xyz[0],nadd,add,nbrks,brks);
  printf(" %4.3f \n",prob);
#endif

#if !DO_NOT_WRITE || 1
  write_ga();
#endif

  ga_avail = 1;

  return;
}

void ZStruct::read_shuttles_data()
{
  if (nshuttle<1)
  {
    printf(" no shuttles, setting wshpair to -1 (niso: %i) \n",niso);
    update_wpair_mem(niso,niso+1);
    for (int i=0;i<niso;i++)
      wshpair[i] = -1;

    return;
  }

  string shfile = "GSMDATA_SH";

  ifstream output(shfile.c_str(),ios::in);
  if (!output)
  {
    printf(" couldn't find %s file \n",shfile.c_str());
    for (int i=0;i<niso;i++)
      wshpair[i] = -1;
    return;
  }

  string line;
  //vector<string> tok_line;

  int nlines = 0;
  while(!output.eof())
  {
    getline(output,line);
    nlines++;
  }
  output.close();

  printf("  found %i lines of shuttles data \n",nlines-1);

  update_wpair_mem(niso,nlines+1);
  ifstream output2(shfile.c_str(),ios::in);
  if (!output2)
  {
    printf(" couldn't find %s file \n",shfile.c_str());
    return;
  }

  int n = 0;
  while(!output2.eof())
  {
    getline(output2,line);
    //tok_line = StringTools::tokenize(line, " \t");
    //printf(" value: %i \n",atoi(tok_line[0].c_str()));
    //wshpair[n++] = atoi(tok_line[0].c_str());
    wshpair[n++] = atoi(line.c_str());
  }
  output2.close();

#if 0
  printf("  shuttles over nlines:");
  for (int i=0;i<nlines-1;i++)
    printf(" %2i",wshpair[i]);
  printf("\n");
#endif

  return;
}


int ZStruct::read_gsm_data()
{
  printf("\n in read_gsm_data \n");

  read_shuttles_data();

  string gsmfile = "GSMDATA";
  ifstream output(gsmfile.c_str(),ios::in);
  if (!output)
  {
    printf(" couldn't find %s file \n",gsmfile.c_str());
    return 0;
  }

  int nfound = 0;
  int nfoundtotal = 0;
  int npairp = npair;
  int wp = -1;
  int wr = -1;
  int ws1, ws2;
  int start = -1;
  int end = 0;
  int isPair = 0;
  string line;
  vector<string> tok_line;
  while(!output.eof())
  {
    getline(output,line);
    //cout << "  RR " << line << endl;

    if (line.find("REACT")!=string::npos)
    { 
      tok_line = StringTools::tokenize(line, " \t");
      if (tok_line.size()<4) 
      {
        printf(" ERROR: REACT line should have 3 entries \n");
        printf(" %s \n",line.c_str());
        exit(1);
      }
      ws1   = atoi(tok_line[1].c_str());
      ws2   = ws1;
      start = atoi(tok_line[2].c_str());
      end   = atoi(tok_line[3].c_str());
      if (end!=-1) end++;
      isPair = 0;
    }

    if (line.find("PAIR")!=string::npos)
    { 
      tok_line = StringTools::tokenize(line, " \t");
      if (tok_line.size()<5) 
      {
        printf(" ERROR: PAIR line should have 4 entries \n");
        printf(" %s \n",line.c_str());
        exit(1);
      }
      ws1   = atoi(tok_line[1].c_str());
      ws2   = atoi(tok_line[2].c_str());
      start = atoi(tok_line[3].c_str());
      end   = atoi(tok_line[4].c_str());
      if (end!=-1) end++;
      isPair = 1;
    }
    if ((start>-1 && end>0) || end==-1)
    {
      printf("\n  found start: %i end: %i for reactant pair: %i %i isPair: %i \n",start,end,ws1,ws2,isPair);

      string nstr = StringTools::int2str(wp+1,4,"0");
      string filename = "scratch/rpair"+nstr+".xyz";
      if (!isPair) 
      {
        nstr = StringTools::int2str(ws1+1,1,"0");
        filename = "react"+nstr+".xyz";
      }

      struct stat sts;
      if (stat(filename.c_str(), &sts) != -1)
      {
        int store = 1;
        if (isPair)
        {
          wp++;

          update_pairs_mem(wp+1,npairp);
 
         //assign_sequential_pairs(wp);
          pairs[2*wp+0] = ws1;
          pairs[2*wp+1] = ws2;
          pairsnum[2*wp+0] = start;
          pairsnum[2*wp+1] = end;
          if (end==-1) pairsnum[2*wp+1] = start;
          npair = wp+1;
          fill_elistrp();
        }
        else
        { 
          wr++;
          update_react_mem();
          reactnum[2*wr+0] = start;
          reactnum[2*wr+1] = end;
          if (end==-1) reactnum[2*wr+1] = start;
          reactrun[wr] = 1;
        }

        read_gsm_files(start,end,NULL);

        if (end!=-1)
          update_wpair_mem(start,end);
        if (isPair)
        {
          for (int i=start;i<end;i++)
            wpair[i] = wp;
          get_elistcmin(wp,end);
        }
        else
        {
          for (int i=start;i<end;i++)
            wpair[i] = -wr-1;
        }

        int wpr = wp;
        if (!isPair) wpr = -wr-1;
        gsm_analyze(wpr,start,end,filename);

        nfound = end - start; 
        if (end==-1) nfound = 0;
        printf("   nfound: %i \n",nfound);
        nfoundtotal += nfound;
      }
      else
        printf(" couldn't find %s \n",filename.c_str());

      start = -1;
      end = 0;
    }
  } //while !eof

  output.close();

#if 1
  printf("\n  after read, pairs: %i \n",npair);
  for (int i=0;i<npair;i++)
    printf("   pair (%3i): %2i %2i range: %4i %4i \n",i,pairs[2*i+0],pairs[2*i+1],pairsnum[2*i+0],pairsnum[2*i+1]);
  printf("\n");
#endif

  printf(" done with reading gsm data: %3i \n",nfoundtotal); fflush(stdout);

#if !DO_NOT_WRITE && 0
  printf("   writing final.xyz \n");
  rxns2.save_structures();
#endif


  return nfoundtotal;
}


void ZStruct::gsm_analyze(int wp, int first, int last, string filename)
{
  printf("\n in gsm_analyze for pair %2i (%3i --> %3i) \n",wp,first,last); fflush(stdout);

  if (last<=first)
  {
    printf("  no analysis possible: %i isomers found \n",last);
    return;
  }
 
 //assuming only one shuttle for now
  int natomssh = 0;
  //for (int i=first;i<last;i++)
  //  printf(" wshpair[%i]: %i \n",i,wshpair[i]);
  for (int i=first;i<last;i++)
  if (wshpair[i]>-1)
  {
    if (nshuttle<=wshpair[i])
    {
      printf("  shuttle %i not found! \n",wshpair[i]);
      exit(1);
    }
    int natoms0 = icshuttle[wshpair[i]].natoms;
    natomssh = natoms0;
    break;
  }
  printf("  file: %s natomssh: %i \n",filename.c_str(),natomssh);

  ICoord ic1;
  ic1.init(filename);

  int natoms1 = ic1.natoms;
  if (natoms1<1)
  {
    printf(" WARNING: no atoms \n");
    return;
  }
  int natoms2 = natoms1 + natomssh;

  ICoord ic2;
  ic2.alloc(natoms2);
  ic2.natoms = natoms2;
  for (int i=0;i<3*natoms1;i++)
    ic2.coords[i] = ic1.coords[i];
  for (int i=0;i<natoms1;i++)
  {
    ic2.anames[i] = ic1.anames[i];
    ic2.anumbers[i] = ic1.anumbers[i];
  }


  int nbrks,nadd;
  int* add = new int[16];
  int* brks = new int[16];
  int maxnodes = 30;
  double** xyz = new double*[maxnodes];
  for (int i=0;i<maxnodes;i++)
    xyz[i] = new double[3*natoms2];
  for (int i=0;i<maxnodes;i++)
  for (int j=0;j<3*natoms2;j++)
    xyz[i][j] = 0.;
  double* energies = new double[maxnodes];
  for (int i=0;i<maxnodes;i++)
    energies[i] = 0.;

  //read in string
  //intended vs actual rxn coordinate
  for (int i=first;i<last;i++)
  if (!rsaved[i])
  {
    //printf("\n\n");
    int success = read_ISOMER(i,nadd,add,nbrks,brks);
    int wsh = wshpair[i];

    int nnodes;
    if (wsh==-1)
      nnodes = read_string(i,energies,natoms1,xyz);
    else
      nnodes = read_string(i,energies,natoms2,xyz);
    if (nnodes<0)
    {
      wsh = wshpair[i] = 0; //shuttle fix
      success = 0;
    }

    if (success && nnodes)
    {
      //printf("  now analyzing and storing %3i \n",i); fflush(stdout);
      if (wsh==-1)
      {
        rsaved[i] = analyze_path(i,wp,nadd,add,nbrks,brks,nnodes,energies,xyz,ic1);
        store_xyz(i,nnodes,natoms1,ic1.anames,ic1.anumbers,xyz);
      }
      else
      {
//#if ISOMER_DB
//        printf("  shouldn't be here in non-shuttle \n");
//        exit(1);
//#endif
        for (int j=0;j<natomssh;j++)
        {
          ic2.anames[natoms1+j] = icshuttle[wsh].anames[j];
          ic2.anumbers[natoms1+j] = icshuttle[wsh].anumbers[j];
        }
        for (int j=0;j<3*natomssh;j++)
          ic2.coords[3*natoms1+j] = icshuttle[wsh].coords[j];
#if ISOMER_DB
        ic2.ic_create();
#endif

        rsaved[i] = analyze_path(i,wp,nadd,add,nbrks,brks,nnodes,energies,xyz,ic2);
        store_xyz(i,nnodes,natoms2,ic2.anames,ic2.anumbers,xyz);
      }
    }
    else
      printf("   insufficient data for structure %3i \n",i);
  }

  fflush(stdout);

  ic1.freemem();
  ic2.freemem();
  delete [] add;
  delete [] brks;
  for (int i=0;i<maxnodes;i++)
    delete [] xyz[i];
  delete [] xyz;
  delete [] energies;

  return;
}


void ZStruct::get_qrc(int ws, double* qrc1)
{
 //gather charge information from reactants
  int wp0 = wpair[ws];
  int wp = wp0;
  if (wp<0) wp = -wp - 1;
  int wr1, wr2;
  if (wp0>-1)
  {
    wr1 = pairs[2*wp+0];
    wr2 = pairs[2*wp+1];
  }
  else
  {
    wr1 = wp;
    wr2 = wp;
  }
  int wsh = wshpair[ws];
  int wshr = 0;
  if (wsh>-1)
    wshr = wshuttle[wsh];
  int naa = natomsr[wr1];
  int nab = natomsr[wr2];
  int nash = natomsr[wshr];
  if (wp0<0) nab = 0;
  if (wsh<0) nash = 0;
  for (int i=0;i<naa+nab+nash;i++)
    qrc1[i] = 9999.;
  for (int i=0;i<naa;i++)
    qrc1[i] = qr[wr1][i];
  for (int i=0;i<nab;i++)
    qrc1[naa+i] = qr[wr2][i];
  if (wsh>-1)
  for (int i=0;i<nash;i++)
    qrc1[naa+nab+i] = qr[wshr][i];

#if 0
  printf("  wr1: %i wr2: %i wsh: %i wshr: %i naa: %i nab: %i nash: %i \n",wr1,wr2,wsh,wshr,naa,nab,nash);
  printf("  get_qrc:");
  for (int i=0;i<naa+nab+nash;i++)
    printf(" %4.2f",qrc1[i]);
  printf("\n");
#endif

  return;
}

void ZStruct::store_xyz(int ws, int nnodes, int natoms0, string* anames0, int* anumbers0, double** xyz0)
{
  //allocates and stores, one allocated at a time
  //pointer array for each ws allocated all at once

  double* xyzr0 = xyz0[0];
  double* xyzp0 = xyz0[nnodes-1];

  int niso1 = niso;
  if (niso<=ws || niso==0)
    niso1 = ws + 1;

  double* qrc0 = new double[natoms0];
  get_qrc(ws,qrc0);

  //printf(" in store_xyz, niso: %i niso1: %i nisoa: %i \n",niso,niso1,nisoa);

  if (niso1>nisoa)
  {
    int* natoms1 = new int[niso1];
    string** anames1 = new string*[niso1];
    int** anumbers1 = new int*[niso1];
    double** xyzr1 = new double*[niso1];
    double** xyzp1 = new double*[niso1];
    //double** xyzt1 = new double*[niso1];
    double** qrc1 = new double*[niso1];
    for (int i=0;i<nisoa;i++)
    {
      natoms1[i] = natoms[i];
      anames1[i] = anames[i];
      anumbers1[i] = anumbers[i];
      xyzr1[i] = xyzr[i];
      xyzp1[i] = xyzp[i];
      //xyzt1[i] = xyzt[i];
      qrc1[i] = qrc[i];
    }
    for (int i=nisoa;i<niso1;i++)
    {
      natoms1[i] = 0;
      anames1[i] = NULL;
      anumbers1[i] = NULL;
      xyzr1[i] = NULL;
      xyzp1[i] = NULL;
      //xyzt1[i] = NULL;
      qrc1[i] = NULL;
    }
    if (nisoa>0)
    {
      delete [] natoms;
      delete [] anames;
      delete [] anumbers;
      delete [] xyzr;
      delete [] xyzp;
      //delete [] xyzt;
      delete [] qrc;
    }

#if 1
    if (anames1[ws]!=NULL) delete [] anames1[ws];
    if (anumbers1[ws]!=NULL) delete [] anumbers1[ws];
    if (xyzr1[ws]!=NULL) delete [] xyzr1[ws];
    if (xyzp1[ws]!=NULL) delete [] xyzp1[ws];
    //if (xyzt1[ws]!=NULL) delete [] xyzt1[ws];
    if (qrc1[ws]!=NULL) delete [] qrc1[ws];
#endif

    natoms1[ws] = 0;
    anames1[ws] = new string[natoms0];
    anumbers1[ws] = new int[natoms0];
    xyzr1[ws] = new double[3*natoms0];
    xyzp1[ws] = new double[3*natoms0];
    //xyzt1[ws] = new double[3*natoms0];
    qrc1[ws] = new double[natoms0];

    nisoa = niso1;
    natoms1[ws] = natoms0;
    natoms = natoms1;
    anames = anames1;
    anumbers = anumbers1;
    xyzr = xyzr1;
    xyzp = xyzp1;
    //xyzt = xyzt1;
    qrc = qrc1;
  } //if allocated more memory

  if (anames[ws]==NULL)
    anames[ws] = new string[natoms0];
  if (anumbers[ws]==NULL)
    anumbers[ws] = new int[natoms0];
  if (xyzr[ws]==NULL)
    xyzr[ws] = new double[3*natoms0];
  if (xyzp[ws]==NULL)
    xyzp[ws] = new double[3*natoms0];
  //if (xyzt[ws]==NULL)
    //xyzt[ws] = new double[3*natoms0];
  if (qrc[ws]==NULL)
    qrc[ws] = new double[natoms0];

  natoms[ws] = natoms0;
  for (int i=0;i<natoms0;i++)
    anames[ws][i] = anames0[i];
  for (int i=0;i<natoms0;i++)
    anumbers[ws][i] = anumbers0[i];
  for (int i=0;i<3*natoms0;i++)
    xyzr[ws][i] = xyzr0[i];
  for (int i=0;i<3*natoms0;i++)
    xyzp[ws][i] = xyzp0[i];
  //for (int i=0;i<3*natoms0;i++)
  //  xyzt[ws][i] = xyzt0[i];
  for (int i=0;i<natoms0;i++)
    qrc[ws][i] = qrc0[i];

  delete [] qrc0;

  return;
}

#if 0
void ZStruct::store_xyz(int ws, int nnodes, int natoms0, string* anames0, int* anumbers0, double** xyz0)
{
 //note this function allocates up to niso, all with natoms0 size..

  double* xyzr0 = xyz0[0];
  double* xyzp0 = xyz0[nnodes-1];

  int niso1 = niso;
  if (niso<ws || niso==0)
    niso1 = ws + 1;

  //printf(" in store_xyz, niso: %i niso1: %i nisoa: %i \n",niso,niso1,nisoa);

  if (niso1>nisoa)
  {
    int* natoms1 = new int[niso1];
    string** anames1 = new string*[niso1];
    int** anumbers1 = new int*[niso1];
    double** xyzr1 = new double*[niso1];
    double** xyzp1 = new double*[niso1];
    //double** xyzt1 = new double*[niso1];
    for (int i=0;i<nisoa;i++)
    {
      natoms1[i] = natoms[i];
      anames1[i] = anames[i];
      anumbers1[i] = anumbers[i];
      xyzr1[i] = xyzr[i];
      xyzp1[i] = xyzp[i];
      //xyzt1[i] = xyzt[i];
    }
    for (int i=nisoa;i<niso;i++)
    {
      natoms1[i] = 0;
      anames1[i] = new string[natoms0];
      anumbers1[i] = new int[natoms0];
      xyzr1[i] = new double[3*natoms0];
      xyzp1[i] = new double[3*natoms0];
      //xyzt1[i] = new double[3*natoms0];
    }
    if (nisoa>0)
    {
      delete [] natoms;
      delete [] anames;
      delete [] anumbers;
      delete [] xyzr;
      delete [] xyzp;
      //delete [] xyzt;
    }
    nisoa = niso1;
    natoms1[ws] = natoms0;
    natoms = natoms1;
    anames = anames1;
    anumbers = anumbers1;
    xyzr = xyzr1;
    xyzp = xyzp1;
    //xyzt = xyzt1;
  }
  for (int i=0;i<natoms0;i++)
    anames[ws][i] = anames0[i];
  for (int i=0;i<natoms0;i++)
    anumbers[ws][i] = anumbers0[i];
  for (int i=0;i<3*natoms0;i++)
    xyzr[ws][i] = xyzr0[i];
  for (int i=0;i<3*natoms0;i++)
    xyzp[ws][i] = xyzp0[i];
  //for (int i=0;i<3*natoms0;i++)
  //  xyzt[ws][i] = xyzt0[i];

  return;
}
#endif


//CPMZ energies array is from stringfile.xyz, using energies from tsq file
int ZStruct::analyze_path(int ws, int wp, int nadd, int* add, int nbrks, int* brks, int nnodes, double* energies, double** xyz, ICoord ic1)
{
  int saved = 1;
  int nmax;
  double emax = get_tse(nnodes,nmax,energies);
  if (emax>300) return 0;
  double* E = new double[3];

  //printf("      analyze_path for %4i \n",ws); fflush(stdout);
  if (nmax<0 || nmax>nnodes)
  {
    printf(" ERROR: problem in nmax: %i \n",nmax);
    return 0;
  }

 //reading energies from stored GSM data
  double esh = 0.;
  int wsh = wshpair[ws];
  if (wsh>-1) esh = elistr[wshuttle[wsh]];
  double emin;
  if (wp>=0)
  {
    if (wsh==-1)
    {
      emin = elistcmin[wp]*627.5;
  //    emin = (elistrp[wp] + esh)*627.5;
    }
    else
    {
      emin = elistcminsh[wp] * 627.5;
  //    emin = (elistrp[wp] + esh)*627.5;
    }
  }
  else 
    emin = elistc[ws]*627.5;

  E[0] = 0.;
  E[1] = elistp[ws]*627.5 - emin;
  E[2] = elistts[ws]*627.5 - emin;

  //if (elistrp!=NULL) printf(" elistcmin[wp]: %6.5f elistrp[wp]: %6.5f esh: %6.5f wshpair[ws]: %i \n",elistcmin[wp],elistrp[wp],esh,wshpair[ws]);
  //printf(" E(%i): %5.1f %5.1f %5.1f  nmax: %2i \n",ws,E[0],E[1],E[2],nmax); fflush(stdout);

  emax = E[2];
  if (E[1]>emax) emax = E[1];

  double* qrc1 = new double[ic1.natoms];
  get_qrc(ws,qrc1);

  if (emax > -0.1 && emax < 250. && !close_val(emax,111.1,0.05))
  {
#if USE_GA || USE_DT || 1
    //int error = rxns2.add_ts(ws,ic1.natoms,ic1.anames,ic1.anumbers,xyz[0],xyz[nnodes-1],xyz[nmax],E,nadd,add,nbrks,brks);
    //printf(" about to add: %i \n",ws); fflush(stdout);
#if !TEST_GA
    int error = rxns2.add_ts_xyz(ws,ic1.natoms,ic1.anames,ic1.anumbers,xyz[0],xyz[nnodes-1],xyz[nmax],E,qrc1);
#else
    int error = rxns2.add_ts_xyz_test(ws,ic1.natoms,ic1.anames,ic1.anumbers,xyz[0],xyz[nnodes-1],xyz[nmax],E,qrc1,nadd,add,nbrks,brks);
#endif
#else
    int error = 0; //CPMZ debug
#endif
    if (error)
    {
      //printf(" didn't add %2i \n",ws);
      saved = 0;
    }

//CPMZ here
#if ISOMER_DB
 //passing string information, initial XYZ
//  ic1.coord_num();
  store_isomers(ws,wp,nadd,nbrks,add,brks,E,ic1.natoms,ic1.anames,ic1.anumbers,xyz[0],qrc1,ic1.coordn);
#endif

  }
  else
  {
    saved = 0;
    //printf(" too high energy or SCF failure, didn't add %2i to db \n",ws);
  }

  delete [] qrc1;
  delete [] E;

  return saved;
}


double ZStruct::get_tse(int nnodes, int& nmax, double* energies)
{
  double emax = -1000.;
  double* E = new double[nnodes];
  for (int i=0;i<nnodes;i++)
    E[i] = energies[i] - energies[0];

  for (int i=0;i<nnodes;i++)
  if (E[i]>emax)
  {
    nmax = i;
    emax = E[i];
  }
  for (int i=0;i<nnodes;i++)
  if (close_val(E[i],111.111,0.005))
  {
    nmax = i;
    emax = 1000.; 
    break;
  }

  delete [] E;

  return emax;
}


void ZStruct::dft_para_q(int niso1)
{
  printf("\n doing dft sp for charges \n"); fflush(stdout);

 //lost memory! (qp is double**)
  if (qp!=NULL) delete [] qp;
  qp = new double*[niso1];

  int natomsmax = 0;
  for (int i=0;i<niso1;i++)
  if (natomsmax<natoms[i])
    natomsmax = natoms[i];

#if !USE_MOPAC
  DFT dft1;
  dft1.alloc(natomsmax);

  printf(" DFT preparing structures for sp ");
  for (int i=0;i<niso1;i++)
  {
    printf(" %i",i);
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/spdft"+nstr;
    dft1.reset(natoms[i],anumbers[i],anames[i],xyzp[i]);
#if !DO_NOT_WRITE
    dft1.opt_dnr(filename);
#endif
  }
  printf("\n");

  int nfound = 0;
  for (int i=0;i<nreact;i++)
  {
    string nstr=StringTools::int2str(i,1,"0"); 
    string dftfile_string = "scratch/dftspdone"+nstr;
    struct stat sts;
    if (stat(dftfile_string.c_str(), &sts) != -1)
    {
      printf(" done with %2i",i); fflush(stdout);
      nfound++;
    }
  }
  printf("\n");

  if (nfound<nreact)
  {
    ofstream cmdfile;
    string cmdfile_string = "scratch/go_dft_sp";

    cmdfile.open(cmdfile_string.c_str());
    printf("\n printing job array numbers for DFT SP \n");
    //cmdfile << "#PBS -t ";
    cmdfile << "#SBATCH --array= ";
    if (niso1-1>0)
      cmdfile << "0-" << niso1-1 << endl;
    else
      cmdfile << "0" << endl;
    cmdfile.close();
    printf("\n");
#if !SKIPDFT
    string cmd = "./qmakedsp";
    system(cmd.c_str());
    cmd = "qsub "+cmdfile_string+".qsh";
    system(cmd.c_str());
#endif

    int qnotdone = true;
    int* dftdone = new int[niso1+1];
    for (int i=0;i<niso1;i++) dftdone[i]=0;

#if !SKIPDFT
    int max_wait = MAX_TIME_WAIT;
    int tc = 0;
    do {
      tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
      sleep(10);
      qnotdone = false;
      for (int i=0;i<nreact;i++)
      {
        if (!dftdone[i])
        {
          string nstr=StringTools::int2str(i,1,"0");
          string dftfile_string = "scratch/dftspdone"+nstr;
          struct stat sts;
          if (stat(dftfile_string.c_str(), &sts) != -1)
          {
            dftdone[i]=1;
            printf(" done with %2i",i); fflush(stdout);
          }
          else
            qnotdone = true;
        } // if not done
      } //loop i over nts
    } while (qnotdone);
    printf("\n");
    sleep(5);
#endif
    delete [] dftdone;
  } //if nfound<nreact
  else
    printf(" not submitting, already done! \n");

  printf(" reading SP data \n"); fflush(stdout);
  for (int i=0;i<niso1;i++)
  {
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/spdft"+nstr;
    dft1.reset(natoms[i],anumbers[i],anames[i],xyzp[i]);
    double ep = dft1.get_opt_energy(filename);
    //dft1.get_structure(filename,icr[i].coords);
//CPMZ  here
    if (qp[i]!=NULL) delete [] qp[i];
    qp[i] = new double[natoms[i]];
    dft1.get_charges(filename,qp[i]);

    printf("  q read in from DFT %i  E is %7.6f \n",i,ep);
    //print_xyz_gen(natoms[i],anames[i],xyzp[i]);
  }

  dft1.freemem();
#else
  printf("\n using MOPAC instead of DFT \n");
  Mopac mopt;
  mopt.alloc(natomsmax);
  for (int i=0;i<nreact;i++)
  {
    printf(" reactant %i \n",i); fflush(stdout);
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/reactmop"+nstr;
//    mopt.reset(icr[i].natoms,icr[i].anumbers,icr[i].anames,icr[i].coords);
    mopt.reset(natoms[i],anumbers[i],anames[i],xyzp[i]); //not sure about p
    double energy = mopt.opt(filename)/627.5;
    printf(" energy: %8.6f \n",energy);
    double ep = energy;

    if (qp[i]!=NULL) delete [] qp[i];
    qp[i] = new double[natoms[i]];
    for (int j=0;j<natoms[i];j++) qp[i][j] = 0.;
  }
#endif

//  for (int i=0;i<niso1;i++)
//    printf(" E[%i]: %1.3f",i,elistr[i]);
//  printf("\n");


  return;
}

void ZStruct::dft_para(int nreact, ICoord* icr)
{
  printf("\n doing dft opt \n"); fflush(stdout);
#if USE_ASE_GAUSSIAN
  dft_para_ase_gaussian(nreact,icr);
#else
  dft_para_qchem(nreact,icr);
#endif

}


void ZStruct::write_charge(int q1, string nstr)
{
  string filename = "scratch/chg"+nstr;
  ofstream chgfile;
  chgfile.open(filename.c_str());

  string chgstr = StringTools::int2str(abs(q1),1,"0");
  if (q1<0) chgfile << "-";
  chgfile << chgstr << endl;

  chgfile.close();

//  printf("  writing %2i as %s to file %s \n",q1,chgstr.c_str(),filename.c_str());

  return;
}

void ZStruct::dft_para_ase_gaussian(int nreact, ICoord* icr)
{
  printf("  using ASE or g09 \n"); fflush(stdout);

  if (qr!=NULL) delete [] qr;
  qr = new double*[nreact];
  for (int i=0;i<nreact;i++)
    qr[i] = NULL;

  int natomsmax = 0;
  for (int i=0;i<nreact;i++)
  if (natomsmax<icr[i].natoms)
    natomsmax = icr[i].natoms;

#if !USE_MOPAC
  //DFT dft1;
  //dft1.alloc(natomsmax);
  int type = USE_ASE_GAUSSIAN;

  printf(" DFT preparing structure for opt ");
  for (int i=0;i<nreact;i++)
  {
    printf(" %i",i);
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/reactdft.xyz"+nstr;
#if !DO_NOT_WRITE
    if (type==2)
      icr[i].print_xyz_save_no_natoms(filename);
    else
      icr[i].print_xyz_save(filename);
    write_charge(icr[i].q1,nstr);
#endif
  }
  printf("\n");
  fflush(stdout);

  int nfound = 0;
  for (int i=0;i<nreact;i++)
  {
    string nstr=StringTools::int2str(i,1,"0"); 
    string dftfile_string = "scratch/dftdone"+nstr;
    struct stat sts;
    if (stat(dftfile_string.c_str(), &sts) != -1)
    {
      printf(" done with %2i",i); fflush(stdout);
      nfound++;
    }
  }
  printf("\n");

  if (nfound<nreact)
  {
    ofstream cmdfile;
    string cmdfile_string = "scratch/go_dft";

    cmdfile.open(cmdfile_string.c_str());
    printf("\n printing job array numbers for DFT OPT \n");
    //cmdfile << "#PBS -t ";
    cmdfile << "#SBATCH --array= ";
    if (nreact-1>0)
      cmdfile << "0-" << nreact-1 << endl;
    else
      cmdfile << "0" << endl;
    cmdfile.close();
    printf("\n");
#if !SKIPDFT
    string cmd = "./qmakeda";
    system(cmd.c_str());
    cmd = "qsub "+cmdfile_string+".qsh";
    system(cmd.c_str());
#endif

    int qnotdone = true;
    int* dftdone = new int[nreact+1];
    for (int i=0;i<nreact;i++) dftdone[i]=0;

#if !SKIPDFT
    int max_wait = MAX_TIME_WAIT;
    int tc = 0;
    do {
      tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
      sleep(10);
      qnotdone = false;
      for (int i=0;i<nreact;i++)
      {
        if (!dftdone[i])
        {
          string nstr=StringTools::int2str(i,1,"0");
          string dftfile_string = "scratch/dftdone"+nstr;
          struct stat sts;
          if (stat(dftfile_string.c_str(), &sts) != -1)
          {
            dftdone[i]=1;
            printf(" done with %2i",i); fflush(stdout);
          }
          else
            qnotdone = true;
        } // if not done
      } //loop i over nts
    } while (qnotdone);
    printf("\n");
    sleep(5);
#endif
    delete [] dftdone;
  } //if nfound<nreact
  else
    printf(" not submitting, already done! \n");

  printf(" reading opt'd geoms \n"); fflush(stdout);
  for (int i=0;i<nreact;i++)
  {
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/reactdft.xyzo"+nstr;
    //read in optimized structure and energy
    if (type==2)
    {
      filename = "scratch/reactdft.gopt"+nstr;
      elistr[i] = read_chk_file(filename,icr[i].natoms,icr[i].coords);
    }
    else
      elistr[i] = read_xyz_file(filename,icr[i].natoms,icr[i].coords);

    if (qr[i]!=NULL) delete [] qr[i];
    qr[i] = new double[icr[i].natoms];
    for (int j=0;j<icr[i].natoms;j++) qr[i][j] = 0.;
//    dft1.get_charges(filename,qr[i]);

    printf("  xyz read in from DFT %i  E is %7.6f \n",i,elistr[i]);
    //print_xyz_gen(icr[i].natoms,icr[i].anames,icr[i].coords);
  }

  //dft1.freemem();
#else
  printf("\n using MOPAC instead of DFT \n");
  Mopac mopt;
  mopt.alloc(natomsmax);
  for (int i=0;i<nreact;i++)
  {
    printf(" reactant %i \n",i); fflush(stdout);
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/reactmop"+nstr;
    mopt.reset(icr[i].natoms,icr[i].anumbers,icr[i].anames,icr[i].coords);
    double energy = mopt.opt(filename)/627.5;
    printf(" energy: %8.6f \n",energy);
    elistr[i] = energy;

    qr[i] = new double[icr[i].natoms];
    for (int j=0;j<icr[i].natoms;j++) qr[i][j] = 0.;
  }
#endif

  for (int i=0;i<nreact;i++)
    printf(" E[%i]: %1.3f",i,elistr[i]);
  printf("\n");


  return;
}


void ZStruct::dft_para_qchem(int nreact, ICoord* icr)
{
#if !USE_MOPAC
  printf(" using Q-Chem \n"); fflush(stdout);
#endif

  if (qr!=NULL) delete [] qr;
  qr = new double*[nreact];

  int natomsmax = 0;
  for (int i=0;i<nreact;i++)
  if (natomsmax<icr[i].natoms)
    natomsmax = icr[i].natoms;

#if !USE_MOPAC
  DFT dft1;
  dft1.alloc(natomsmax);

  printf(" DFT preparing structure for opt ");
  for (int i=0;i<nreact;i++)
  {
    printf(" %i",i);
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/reactdft"+nstr;
    dft1.reset(icr[i].natoms,icr[i].anumbers,icr[i].anames,icr[i].coords);
#if !DO_NOT_WRITE
    dft1.opt_dnr(filename);
#endif
  }
  printf("\n");

  int nfound = 0;
  for (int i=0;i<nreact;i++)
  {
    string nstr=StringTools::int2str(i,1,"0"); 
    string dftfile_string = "scratch/dftdone"+nstr;
    struct stat sts;
    if (stat(dftfile_string.c_str(), &sts) != -1)
    {
      printf(" done with %2i",i); fflush(stdout);
      nfound++;
    }
  }
  printf("\n");

  if (nfound<nreact)
  {
    ofstream cmdfile;
    string cmdfile_string = "scratch/go_dft";

    cmdfile.open(cmdfile_string.c_str());
    printf("\n printing job array numbers for DFT OPT \n");
    //cmdfile << "#PBS -t ";
    cmdfile << "#SBATCH --array= ";
    if (nreact-1>0)
      cmdfile << "0-" << nreact-1 << endl;
    else
      cmdfile << "0" << endl;
    cmdfile.close();
    printf("\n");
#if !SKIPDFT
    string cmd = "./qmaked";
    system(cmd.c_str());
    cmd = "qsub "+cmdfile_string+".qsh";
    system(cmd.c_str());
#endif

    int qnotdone = true;
    int* dftdone = new int[nreact+1];
    for (int i=0;i<nreact;i++) dftdone[i]=0;

#if !SKIPDFT
    int max_wait = MAX_TIME_WAIT;
    int tc = 0;
    do {
      tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
      sleep(10);
      qnotdone = false;
      for (int i=0;i<nreact;i++)
      {
        if (!dftdone[i])
        {
          string nstr=StringTools::int2str(i,1,"0");
          string dftfile_string = "scratch/dftdone"+nstr;
          struct stat sts;
          if (stat(dftfile_string.c_str(), &sts) != -1)
          {
            dftdone[i]=1;
            printf(" done with %2i",i); fflush(stdout);
          }
          else
            qnotdone = true;
        } // if not done
      } //loop i over nts
    } while (qnotdone);
    printf("\n");
    sleep(5);
#endif
    delete [] dftdone;
  } //if nfound<nreact
  else
    printf(" not submitting, already done! \n");

  printf(" reading opt'd geoms \n"); fflush(stdout);
  for (int i=0;i<nreact;i++)
  {
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/reactdft"+nstr;
    dft1.reset(icr[i].natoms,icr[i].anumbers,icr[i].anames,icr[i].coords);
    elistr[i] = dft1.get_opt_energy(filename);
    dft1.get_structure(filename,icr[i].coords);
//CPMZ  here
    qr[i] = new double[icr[i].natoms];
    dft1.get_charges(filename,qr[i]);

    printf("  xyz read in from DFT %i  E is %7.6f \n",i,elistr[i]);
    //print_xyz_gen(icr[i].natoms,icr[i].anames,icr[i].coords);
  }

  dft1.freemem();
#else
  printf("\n using MOPAC instead of DFT \n");
  Mopac mopt;
  mopt.alloc(natomsmax);
  for (int i=0;i<nreact;i++)
  {
    printf(" reactant %i \n",i); fflush(stdout);
    string nstr=StringTools::int2str(i,4,"0"); 
    string filename = "scratch/reactmop"+nstr;
    mopt.reset(icr[i].natoms,icr[i].anumbers,icr[i].anames,icr[i].coords);
    double energy = mopt.opt(filename)/627.5;
    printf(" energy: %8.6f \n",energy);
    elistr[i] = energy;

    qr[i] = new double[icr[i].natoms];
    for (int j=0;j<icr[i].natoms;j++) qr[i][j] = 0.;
  }
#endif

  for (int i=0;i<nreact;i++)
    printf(" E[%i]: %1.3f",i,elistr[i]);
  printf("\n");


  return;
}

void ZStruct::gsm_para(int first, int last)
{
  printf("\n in gsm_para \n"); fflush(stdout);

  ofstream cmdfile;
  string cmdfile_string = "scratch/go_gsm_dft";

  cmdfile.open(cmdfile_string.c_str());
  //cmdfile << "#PBS -t ";
  cmdfile << "#SBATCH --array= ";
  cmdfile << first << "-" << last-1 << endl;
  cmdfile.close();

  printf(" gsm_para range: %i-%i \n",first,last-1); fflush(stdout);
  if (first>last-1) return;

  string cmd = "./qmakegf";
  system(cmd.c_str());
  cmd = "qsub "+cmdfile_string+".qsh";
#if !SKIPDFT && !SKIPGSM
  system(cmd.c_str());
#else
  printf(" not submitting to queue! \n");
#endif

  int qnotdone = true;
  int* gsmdone = new int[last];
  for (int i=0;i<last;i++) gsmdone[i]=1;
  for (int i=first;i<last;i++) gsmdone[i]=0;

#if DO_NOT_WAIT
  int max_wait = 0;
#else
  int max_wait = MAX_TIME_WAIT*10;
#endif

#define READ_ONLY 1

  int ngdone = 0;
  int ngtotal = last - first + 1;
  int tc = 0;
  printf("\n done list:");
  do {
    tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
    sleep(15);
    if (tc%10==0) fflush(stdout);
    qnotdone = false;
    for (int i=first;i<last;i++)
    {
      if (!gsmdone[i])
      {
#if READ_ONLY
        ngdone++;
        gsmdone[i] = 1;
#else
        string nstr=StringTools::int2str(i,1,"0");
        string gsmfile_string = "scratch/gsmdone"+nstr;
        struct stat sts;
        if (stat(gsmfile_string.c_str(), &sts) != -1)
        {
          ngdone++;
          gsmdone[i]=1;
          printf(" %i",i);
        }
        else
          qnotdone = true;
#endif
      } // if not done
    } //loop i over gsm runs
    double ngdoned = ngdone;
    double fdone = ngdoned/ngtotal;
    if (fdone>GSM_FRAC_DONE)
      qnotdone = false;
    if (tc%100==0) printf(" f: %3.2f",fdone);

  } while (qnotdone);
  printf("\n");

  read_gsm_files(first,last,gsmdone);
  
  delete [] gsmdone;

  return;
}


void ZStruct::read_gsm_files(int first, int last, int* gsmdone)
{
  printf("\n  in read_gsm_files: %3i to %3i \n",first,last); fflush(stdout);

  if (last<=first)
  {
    printf("   nothing to read (first: %i last: %i) \n",first,last);
    return;
  }

  update_elist_mem(first, last);

  int delgsmdone = 0; 
  if (gsmdone==NULL)
  {
    gsmdone = new int[last];
    for (int i=0;i<last;i++)
      gsmdone[i] = 1;
    delgsmdone = 1;
  }

  double* energies = new double[30];

  for (int i=first;i<last;i++)
  if (gsmdone[i])
  {
    string nstr = StringTools::int2str(i,4,"0");
    string oname = "scratch/tsq"+nstr+".xyz";
    ifstream output(oname.c_str(),ios::in);
    if (!output)
    {
      oname = "scratch/savetsq/tsq"+nstr+".xyz";
      output.open(oname.c_str());
    }
    string line;
    int foundtsq = 0;
    int natoms1 = 0;
    if (output)
    {
      getline(output,line);
      natoms1 = atoi(line.c_str());

      getline(output,line);
      //cout << " RR " << line << endl;

      vector <string> tok_line;
      tok_line = StringTools::tokenize(line, " \t");

      double energyts = 0.;
      if (tok_line.size()>2)
      {
        energyts = atof(tok_line[0].c_str());
        double energyr = atof(tok_line[1].c_str());
        double energyp = atof(tok_line[2].c_str());
        elistc[i] = energyr;
        elistts[i] = energyts;
        elistp[i] = energyp;
        //printf("  energies c/ts/p (%i): %6.5f %6.5f %6.5f \n",i,energyr,energyts,energyp);
        //energyts = (energyts-energyr)*627.5;
        foundtsq = 1;
      }
      else printf("  no r/p data available for %2i \n",i);
      if (abs(energyts)<0.00001)
      {
        printf("   GSM failed for %4i \n",i);
        energyts = 25000.;
      }
//      else
//        printf(" final GSM energy for %i is: %8.5f \n",i,energyts);

    }//if output exists 
    else
    {
      //printf(" couldn't find %s \n",oname.c_str());
      elistc[i] = 0.;
      elistts[i] = 10000./627.5;
      elistp[i] = 10000./627.5;
    }
    output.close();

    if (foundtsq)
    {
      int nnodes = read_string(i,energies,natoms1,NULL);
#if 0
      printf("   E:");
      for (int j=0;j<nnodes;j++)
        printf(" %4.1f",energies[j]);
      printf("\n");
#endif
      int founduphill = 0;
      double emin = 0.;
      for (int j=1;j<nnodes;j++)
      {
        if (emin>energies[j])
          emin = energies[j];
        if (energies[j]>energies[j-1]+0.2) //going back uphill
        {
          founduphill = 1;
          break;
        }
      }
      if (emin < -15.) emin = -15.; //in case reaction occurs while going downhill
      //printf("   adjusting elistc for %5i by %4.1f \n",i,emin);
      if (founduphill) 
        elistc[i] += emin/627.5; //emin in kcal/mol
    }

  } //loop over read tsq

  delete [] energies;
  if (delgsmdone) delete [] gsmdone;

  return;
}

void ZStruct::get_elistref()
{
  printf("  in get_elistref. niso: %3i npair: %3i \n",niso,npair); fflush(stdout);

  for (int i=0;i<niso;i++)
  if (wpair[i]<0)
  {
    //int wp = wpair[i];
    elistref[i] = elistc[i];
    int refd = formula_match(0,i);
    if (refd && elistc[0]<elistc[i])
      elistref[i] = elistc[0];
  }
  else
    elistref[i] = elistc[i];

  return;
}

int ZStruct::formula_match(int ws1, int ws2)
{
 //does not rotate atom indices
  if (natoms[ws1]!=natoms[ws2])
    return 0;

  int match = 1;
  int natoms1 = natoms[ws1];
  for (int i=0;i<natoms1;i++)
  if (anames[ws1][i]!=anames[ws2][i])
  {
    match = 0;
    break;
  }

  return match;
}


void ZStruct::get_elistcmin(int wp, int niso1)
{
  int pstart = 0;
  int pend = 0;

  get_pair_indices(wp,niso1,pstart,pend);
  printf("   get_elistcmin for wp: %2i, pstart: %i pend: %i \n",wp,pstart,pend);

  int ws = -1;
  double emin = 100.;
  for (int i=pstart;i<pend;i++)
  if (elistc[i]<emin && wshpair[i]==-1)
  {
    emin = elistc[i];
    ws = i;
  }
  elistcmin[wp] = emin;

  int wsh = -1;
  double eminsh = 100.;
  for (int i=pstart;i<pend;i++)
  if (elistc[i]<eminsh && wshpair[i]!=-1)
  {
    eminsh = elistc[i];
    wsh = i;
  }
  elistcminsh[wp] = eminsh;

  printf("    found emin: %8.5f eminsh: %8.5f (structure: %5i) \n",emin,eminsh,ws);

  return;
}


void ZStruct::get_pair_indices(int wp, int niso1, int& pstart, int& pend)
{
//assumes all pair isos are sequential 

  for (int i=0;i<niso1;i++)
  if (wpair[i]==wp)
  {
    pstart = i;
    break;
  }

  for (int i=0;i<niso1;i++)
  if (wpair[i]==wp)
    pend = i;

  return;
}

void ZStruct::update_react_mem()
{
  int* reactnum1 = new int[2*nreact];
  int* reactrun1 = new int[nreact];
  int* natomsr1 = new int[nreact];
  double** qr1 = new double*[nreact];
  double* elistr1 = new double[nreact];
  int* pair_react1 = new int[nreact];

  if (reactnum!=NULL)
  {
    for (int i=0;i<2*nreacta;i++)
      reactnum1[i] = reactnum[i];
  }
  if (reactrun!=NULL)
  {
    for (int i=0;i<nreacta;i++)
      reactrun1[i] = reactrun[i];
    delete [] reactrun;
  }
  if (natomsr!=NULL)
  {
    for (int i=0;i<nreacta;i++)
      natomsr1[i] = natomsr[i];
    delete [] natomsr;
  }
  if (qr!=NULL)
  {
    for (int i=0;i<nreacta;i++)
      qr1[i] = qr[i];
    delete [] qr;
  }
  if (elistr!=NULL)
  {
    for (int i=0;i<nreacta;i++)
      elistr1[i] = elistr[i];
    delete [] elistr;
  }
  if (pair_react!=NULL)
  {
    for (int i=0;i<nreacta;i++)
      pair_react1[i] = pair_react[i];
    delete [] pair_react;
  }

  for (int i=2*nreacta;i<2*nreact;i++)
    reactnum1[i] = -1;
  for (int i=nreacta;i<nreact;i++)
    reactrun1[i] = 0;
  for (int i=nreacta;i<nreact;i++)
    natomsr1[i] = 0;
  for (int i=nreacta;i<nreact;i++)
    qr1[i] = NULL;
  for (int i=nreacta;i<nreact;i++)
    elistr1[i] = 0.;
  for (int i=nreacta;i<nreact;i++)
    pair_react1[i] = 1;

  reactnum = reactnum1;
  reactrun = reactrun1;
  natomsr = natomsr1;
  qr = qr1;
  elistr = elistr1;
  pair_react = pair_react1;

  nreacta = nreact;

  return;
}

void ZStruct::update_pairs_mem(int npairs0, int npairp)
{
  printf("   in update_pairs_mem: %i npairp: %i \n",npairs0,npairp);

  if (npairs0<=npairp) return;

  int* pairs1 = new int[3*npairs0];
  int* pairsnum1 = new int[2*npairs0];
  double* elistcmin1 = new double[npairs0];
  double* elistcminsh1 = new double[npairs0];

  if (pairs!=NULL)
  for (int i=0;i<2*npair;i++)
    pairs1[i] = pairs[i];

  if (pairsnum!=NULL)
  for (int i=0;i<2*npair;i++)
    pairsnum1[i] = pairsnum[i];

  if (elistcmin!=NULL)
  for (int i=0;i<npair;i++)
    elistcmin1[i] = elistcmin[i];

  if (elistcminsh!=NULL)
  for (int i=0;i<npair;i++)
    elistcminsh1[i] = elistcminsh[i];

  if (pairs!=NULL) delete [] pairs;
  if (pairsnum!=NULL) delete [] pairsnum;
  if (elistcmin!=NULL) delete [] elistcmin;
  if (elistcminsh!=NULL) delete [] elistcminsh;

  pairs = pairs1;
  pairsnum = pairsnum1;
  elistcmin = elistcmin1;
  elistcminsh = elistcminsh1;

#if 0
  printf("   after update_pairs_mem \n");
  for (int i=0;i<npair;i++)
    printf("   pair (%i): %i %i range: %i %i \n",i,pairs[2*i+0],pairs[2*i+1],pairsnum[2*i+0],pairsnum[2*i+1]);
#endif

  return;
}

void ZStruct::update_wpair_mem(int nisop, int nisof)
{
  if (nisop+nisof<=wpaira) return;

  printf("   in update_wpair_mem, nisop: %i nisof: %i wpaira: %i \n",nisop,nisof,wpaira); fflush(stdout);

  int* wpair1 = new int[nisop+nisof+1];
  int* wshpair1 = new int[nisop+nisof+1];
  int* rsaved1 = new int[nisop+nisof+1];

  if (wpair!=NULL)
  for (int i=0;i<wpaira;i++)
    wpair1[i] = wpair[i];
  if (wshpair!=NULL)
  for (int i=0;i<wpaira;i++)
    wshpair1[i] = wshpair[i];
  for (int i=wpaira;i<nisop+nisof;i++)
    wshpair1[i] = -1;
  if (rsaved!=NULL)
  for (int i=0;i<wpaira;i++)
    rsaved1[i] = rsaved[i];
  for (int i=wpaira;i<nisop+nisof+1;i++)
    rsaved1[i] = 0;

  if (wpair!=NULL) delete [] wpair;
  if (wshpair!=NULL) delete [] wshpair;
  if (rsaved!=NULL) delete [] rsaved;

  wpair = wpair1;
  wshpair = wshpair1;
  rsaved = rsaved1;

  wpaira = nisop+nisof;
  //printf(" end update_wpair_mem \n"); fflush(stdout);

  return;
}

void ZStruct::update_elist_mem(int start, int end)
{
  printf("   in update_elist_mem: %i %i elista: %i \n",start,end,elista); fflush(stdout);
  if (end<=elista) return;

  double* elistc1 = new double[end+2];
  double* elistts1 = new double[end+2];
  double* elistp1 = new double[end+2];
  double* elistref1 = new double[end+2];

  for (int i=start;i<end;i++) elistc1[i] = 99.;
  if (elistc!=NULL)
  for (int i=0;i<start;i++)
  {
    elistc1[i] = elistc[i];
    elistts1[i] = elistts[i];
    elistp1[i] = elistp[i];
    elistref1[i] = elistref[i];
  }

  if (elistc!=NULL) delete [] elistc;
  if (elistts!=NULL) delete [] elistts;
  if (elistp!=NULL) delete [] elistp;
  if (elistref!=NULL) delete [] elistref;

  elistc = elistc1;
  elistts = elistts1;
  elistp = elistp1;
  elistref = elistref1;

  elista = end+1;

//  elistc[end] = NULL;
//  elistts[end] = NULL;
//  elistp[end] = NULL;
//  elistref[end] = NULL;

  //printf(" end update_elist_mem \n"); fflush(stdout);
  return;
}


int ZStruct::create_isos(ICoord ic1, int nr, int niso1, int* active1, int wsh)
{
  printf("\n in create_isos (wsh: %i) \n",wsh);

  int nfound = 0;

  //for loops to get all possible
  int maxa = MAX_ADD;
  int maxb = MAX_BRK;

  int natoms1 = ic1.natoms;

 //probably don't need this.. instead operate on ic1
  ICoord ic1tm;
  ic1tm.alloc(ic1.natoms);
  ic1tm.reset(ic1.natoms,ic1.anames,ic1.anumbers,ic1.coords);
  ic1tm.q1 = ic1.q1;
  ic1tm.ic_create();
  ic1tm.update_ic();
  int ntmbonds = ic1tm.ic_create_tm();
  for (int i=0;i<ic1.natoms;i++)
  if (ic1.isTM(i))
    ic1.coordn[i] = ic1tm.coordn[i];

#if 1
  if (ntmbonds>0)
  {
    printf("  TM coordn after ic_create_tm: \n");
    for (int i=0;i<ic1.natoms;i++)
    if (ic1.isTM(i))
      printf("   %2s: %2i \n",ic1.anames[i].c_str(),ic1.coordn[i]);
    printf("\n");
  }
#endif

  int nact1 = 0;
  for (int i=0;i<natoms1;i++)
  if (active1[i]) 
    nact1++;
 
  if (nact1>35)
  {
    printf("\n\n WARNING: too many unfrozen atoms! \n");
    printf("   this limit hardcoded: %i max \n",35);
    exit(1);
  }


 //# of possible break combinations
  int* sizeb = new int[maxb];
  sizeb[0] = ic1.nbonds;
  for (int i=1;i<maxb;i++)
    sizeb[i] = sizeb[i-1]*(ic1.nbonds-i);
  for (int i=1;i<maxb;i++)
    sizeb[i] = sizeb[i]/fact(i+1);
  int sizebsum = 0;
  for (int i=0;i<maxb;i++)
    sizebsum += sizeb[i];

  printf("  # of break combinations:");
  for (int i=0;i<maxb;i++)
    printf(" %2i",sizeb[i]);
  printf(" total: %3i \n",sizebsum);

  int* brks = new int[2*sizebsum*maxb];

  //get all breaks, get coordn
  //for all breaks, get adds
  for (int nbrks=0;nbrks<=maxb;nbrks++)
  {
    int nbf = get_breaks(nbrks,brks,ic1,active1);
    for (int nadd=0;nadd<=maxa;nadd++)
    {
      int nif = get_combo_h1(nadd,nbf,nbrks,brks,ic1,ic1tm,nr,active1,wsh);
      nfound += nif;
    }
  }

#if 0
 //need to modify
  if (nreact==1)
  for (int i=0;i<nfound;i++)
    write_initial_xyz(niso1+i,natoms1,ic1.anames,ic1.coords,0);
#endif

  ic1tm.freemem();
  delete [] brks; 


  return nfound;
}



void ZStruct::create_frag_bonds(ICoord ic1)
{
  int nfrags = ic1.nfrags;
  nfragb = 0;

  if (nfrags<2) return;

  printf(" connecting %i fragments \n",nfrags); 

  if (fragb!=NULL) delete [] fragb;
  fragb = new int[nfrags*nfrags];
  for (int i=0;i<nfrags*nfrags;i++) fragb[i] = -1;

  int natoms1 = ic1.natoms;
  for (int f1=0;f1<nfrags;f1++)
  for (int f2=0;f2<f1;f2++)
  {
    int found = 0;
    for (int i=0;i<natoms1;i++) if (!found)
    for (int j=0;j<i;j++) if (!found)
    {
      if (ic1.frags[i]==f1 && ic1.frags[j]==f2)
      {
        printf(" frag bond: %i %i \n",i+1,j+1);
        fragb[2*nfragb+0] = i;
        fragb[2*nfragb+1] = j;
        nfragb++;
        found = 1;
      }
    } //loop i,j over atom pairs
  }

  return;
}


//CPMZ produces duplicates!
 //behavior: only does H transfer to non-H
 //H must be singly attached (check coordn==1 code)
int ZStruct::get_h_transfer(ICoord ic1, int nadd, int* add, int nbrks, int* brksp, int natoms, int* coordn, int* moving, int wsh)
{
  //note: this function designed/tested to work with 2+/1- only
  if (MAX_CHG>=4)
    return nbrks;
//  return nbrks; //debug

  int nbrksp = nbrks;

  //printf(" in get_h_transfer for nadd: %i nbrks: %i natoms: %i \n",nadd,nbrks,natoms);

  //first gather H that add (to heavy atom) but not break
  int nt1 = 0;
  int* t1 = new int[2*nadd];
  int nt2 = 0;
  int* t2 = new int[2*nadd];
  for (int i=0;i<nadd;i++)
  {
    int a1 = add[2*i+0];
    int a2 = add[2*i+1];
    if (ic1.anumbers[a1]==1 && ic1.coordn[a1]==1 && ic1.anumbers[a2]!=1)
    {
      int found = 0;
      for (int j=0;j<2*nbrks;j++)
      if (a1==brksp[j])
      {
        found = 1;
        break;
      }

      if (!found)
      {
        for (int j=0;j<nt1;j++)
        if (a1==t1[j])
          found = 1;
        if (!found)
          t1[nt1++] = a1;
      } // add only H
      else
      {
        found = 0;
        for (int j=0;j<nt2;j++)
        if (a1==t2[j])
          found = 1;
        if (!found)
          t2[nt2++] = a1;
      } //if add/brk H

    } //if H on add list

    if (ic1.anumbers[a2]==1 && ic1.coordn[a2]==1 && ic1.anumbers[a1]!=1)
    {
      int found = 0;
      for (int j=0;j<2*nbrks;j++)
      if (a2==brksp[j])
      {
        found = 1;
        break;
      }

      if (!found)
      {
        for (int j=0;j<nt1;j++)
        if (a2==t1[j])
          found = 1;
        if (!found)
          t1[nt1++] = a2;
      } // add only H
      else
      {
        found = 0;
        for (int j=0;j<nt2;j++)
        if (a2==t2[j])
          found = 1;
        if (!found)
          t2[nt2++] = a2;
      } //if add/brk H

    } //if H on add list

  } //loop i over nadd

 //handle duplicates
  if (nt2>0 && nt1>0 && wsh==-1) //will shuttle first H transfer
  {
    if (t1[0]>t2[0])
      nt1 = 0;
    //printf(" transfers: %i (%s) and %i (%s) \n",t1[0],ic1.anames[t1[0]].c_str(),t2[0],ic1.anames[t2[0]].c_str());
  }

  if (nt1==1)
  for (int i=0;i<nt1;i++)
  {
    int b1 = t1[i];
    int b2 = -1;
    if (ic1.coordn[b1]==1)
    for (int i=0;i<ic1.nbonds;i++)
    {
      if (ic1.bonds[i][0]==b1)
      {
        b2 = ic1.bonds[i][1];
        break;
      }
      else if (ic1.bonds[i][1]==b1)
      {
        b2 = ic1.bonds[i][0];
        break;
      }
    }
    brksp[2*nbrksp+0] = b1;
    brksp[2*nbrksp+1] = b2;
    nbrksp++;
    coordn[b1]--;
    coordn[b2]--;
    moving[b2] = 1;
  }

#if 0
  for (int i=0;i<nadd;i++)
    printf("  add: %i %i  (%s %s) \n",add[2*i+0]+1,add[2*i+1]+1,ic1.anames[add[2*i+0]].c_str(),ic1.anames[add[2*i+1]].c_str());
  for (int i=0;i<nbrksp;i++)
    printf("  brk: %i %i  (%s %s) \n",brksp[2*i+0]+1,brksp[2*i+1]+1,ic1.anames[brksp[2*i+0]].c_str(),ic1.anames[brksp[2*i+1]].c_str());
  for (int i=0;i<nt;i++)
    printf("  atom %i (%s) on add not break list \n",t[i]+1,ic1.anames[t[i]].c_str());
#endif

  delete [] t1;
  delete [] t2;

  return nbrksp;
}



int ZStruct::get_combo_h1(int nadd, int nbf, int nbrks, int* brks, ICoord ic1, ICoord ic1tm, int nr, int* active1, int wsh)
{
  if (nadd+nbf<1) return 0;
  if (!rxns2.tree_screen_0(nadd,nbrks))
  {
    printf("  %i/%i combo screened by tree! \n",nadd,nbrks);
    return 0;
  }

  printf("  in get_combo_h1: %i/%i \n",nadd,nbrks);

  int maxchg = MAX_CHG;
  int minchg = MIN_CHG;
  if (nadd+nbrks > maxchg)
  {
    printf("   connectivity change limit exceeded (%i/%i) \n",nadd+nbrks,maxchg);
    return 0;
  }
  if (nadd+nbrks<minchg)
  {
    printf("   connectivity change too small (%i/%i) \n",nadd+nbrks,minchg);
    return 0;
  }

  //printf(" getting combo with %i atoms \n",ic1.natoms); fflush(stdout);
  //ic1.print_ic();

  int natoms1 = ic1.natoms;
  int* coordn = new int[natoms1];
  int* coordnb = new int[natoms1]; //after breaks are set

  int nfound = 0;
  int* adds = new int[nadd*2+2];

  int* moving = new int[natoms1];
  int* movingb = new int[natoms1];
  nhh = 0;

//  int hone = 0;
//  if (climit_h[1]<2) hone = 1;
  int* brksp = new int[2*(nbrks+2)];

  int* satoms = new int[3];

  int* nangles = new int[10];
  int** angles = new int*[60];
  for (int i=0;i<60;i++)
    angles[i] = new int[3];
  double* anglev = new double[60];
  int* ntor = new int[10];
  int** tor = new int*[40];
  for (int i=0;i<40;i++)
    tor[i] = new int[4];
  double* torv = new double[30];

#if 1
  for (int i=0;i<nbf;i++)
  {
    for (int j=0;j<natoms1;j++) coordnb[j] = ic1.coordn[j];
    for (int j=0;j<natoms1;j++) movingb[j] = 0;

    for (int j=0;j<nbrks;j++)
    {
      int b1 = brks[2*nbrks*i+2*j+0];
      int b2 = brks[2*nbrks*i+2*j+1];
      coordnb[b1]--;
      coordnb[b2]--;
      movingb[b1] = 1;
      movingb[b2] = 1;
      //printf(" b1: %2i b2: %2i \n",b1+1,b2+1);
    }

    if (nadd==0 && nr==1 && wsh==-1)
    {
      int ool = 0;
      for (int j=0;j<natoms1;j++)
      if (active1[j]!=1 && movingb[j])
        ool = 1;
      if (!ool)
      for (int j=0;j<natoms1;j++)
      if (active1[j]==1 && movingb[j])
      {
        //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
        if (coordnb[j]<climit_l[ic1.anumbers[j]] || coordnb[j]>climit_h[ic1.anumbers[j]])
          ool = 1;
      }
      //Note: nothing to handle for TM if only breaking
      if (!ool)
      {
        int ga_pass = ml_eval(ic1,0,NULL,nbrks,&brks[2*nbrks*i]);
        if (ga_pass)
        {
          write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
          write_ISOMER(niso++,0,NULL,nbrks,&brks[2*nbrks*i]);
          nfound++;
        }
      } //if !ool
    } //if nadd==0 && nr==1

    if (nadd==1)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    {
      int ool = 0;
      for (int j=0;j<natoms1;j++) coordn[j] = coordnb[j];
      for (int j=0;j<natoms1;j++) moving[j] = movingb[j];

      if (!ic1.bond_exists(a1,a2))
      {
        coordn[a1]++;
        coordn[a2]++;
        adds[0] = a1;
        adds[1] = a2;
        moving[a1] = 1;
        moving[a2] = 1;
        //printf(" a1: %2i a2: %2i \n",a1+1,a2+1);

        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn);
        if (!ool && nr>1)
          ool = two_frags(a1,a2,ic1);

       //determine if angle driving
        int tm_in = 0;
        if (ic1.isTM(a1)) tm_in++; if (ic1.isTM(a2)) tm_in++;
        for (int j=0;j<2*nbrks;j++) if (ic1.isTM(brks[2*nbrks*i+j])) tm_in++;
        if (!ool && tm_in)
        {
          nfound += tm_angle_drive(nr,nadd,adds,nbrks,&brks[2*nbrks*i],ic1tm,nangles,angles,anglev,ntor,tor,torv,active1);
        }
        else if (!ool)
        {
          if (wsh==-1)
          {
            int ga_pass = ml_eval(ic1,nadd,adds,nbrks,&brks[2*nbrks*i]);
            if (ga_pass)
            {
              if (nr>1)
              {
                align1.add_align(nadd,adds);
                write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
              }
              else
                write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
              write_ISOMER(niso++,nadd,adds,nbrks,&brks[2*nbrks*i]);
              nfound++;
            }
          }
          else if (wsh>-1)
          {
            if (nr>1)
              align1.add_align(nadd,adds);
            int wsh1 = is_shuttle(ic1,nadd,adds,nbrks,&brks[2*nbrks*i],satoms);
            if (wsh==wsh1)
            {
              nfound += get_shuttle(ic1,nr,wsh,nadd,adds,nbrks,&brks[2*nbrks*i],satoms);
            }
          }

        } //if !ool
      } //if a1-a2 bond exists
    } //if nadd1==1, loop over a1,a2


   //CPMZ check me
    if (nadd==2)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    for (int a3=0;a3<=a1;a3++)
    if (active1[a3]==1)
    for (int a4=0;a4<a3;a4++)
    if (active1[a4]==1)
    if (!(a1==a3 && a4>=a2))
    {
      int ool = 0;
      for (int j=0;j<natoms1;j++) coordn[j] = coordnb[j];
      for (int j=0;j<natoms1;j++) moving[j] = movingb[j];

      if (!ic1.bond_exists(a1,a2) && !ic1.bond_exists(a3,a4))
      {
        coordn[a1]++;
        coordn[a2]++;
        coordn[a3]++;
        coordn[a4]++;
        adds[0] = a1;
        adds[1] = a2;
        adds[2] = a3;
        adds[3] = a4;
        moving[a1] = 1;
        moving[a2] = 1;
        moving[a3] = 1;
        moving[a4] = 1;
        //printf("     a1: %2i a2: %2i a3: %2i a4 %2i \n",a1+1,a2+1,a3+1,a4+1);

        if (!ool && nr>1)
          ool = two_frags(a1,a2,a3,a4,ic1);

        int nbrksp = nbrks;
        for (int j=0;j<2*nbrks;j++)
          brksp[j] = brks[2*nbrks*i+j];
#if USE_H_TRANSFER
        if (!ool)
          nbrksp = get_h_transfer(ic1,nadd,adds,nbrks,brksp,natoms1,coordn,moving,wsh);
#endif
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn) + hh_bond(a3,a4,ic1.anumbers,coordn);

        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf("     coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }

       //CPMZ here
       //determine if angle driving
        //printf(" pairs: %i-%i %i-%i TM? %i %i %i %i ",a1,a2,a3,a4,ic1.isTM(a1),ic1.isTM(a2),ic1.isTM(a3),ic1.isTM(a4));
        int tm_in = 0;
        if (ic1.isTM(a1)) tm_in++; if (ic1.isTM(a2)) tm_in++;
        if (ic1.isTM(a3)) tm_in++; if (ic1.isTM(a4)) tm_in++;
        for (int j=0;j<2*nbrks;j++) if (ic1.isTM(brks[2*nbrks*i+j])) tm_in++;
        if (!ool && tm_in)
        {
         //note tm_angle_drive doesn't use brksp (from get_h_transfer)
          nfound += tm_angle_drive(nr,nadd,adds,nbrks,&brks[2*nbrks*i],ic1tm,nangles,angles,anglev,ntor,tor,torv,active1);
        }
        else if (!ool)
        {
          if (wsh==-1)
          {
            int ga_pass = ml_eval(ic1,nadd,adds,nbrksp,brksp);
            if (ga_pass)
            {
              if (nr>1)
              {
                align1.add_align(nadd,adds);
                write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
              }
              else
                write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
              write_ISOMER(niso++,nadd,adds,nbrksp,brksp);
              nfound++;
            }
          } //if !shuttle
          else if (wsh>-1)
          {
            //printf("   not ool, attempting shuttle \n");
            if (nr>1)
              align1.add_align(nadd,adds);
            int wsh1 = is_shuttle(ic1,nadd,adds,nbrksp,brksp,satoms);
            if (wsh==wsh1)
            {
              //printf("   is shuttle true \n");
              nfound += get_shuttle(ic1,nr,wsh,nadd,adds,nbrksp,brksp,satoms);
            }
          }

        } //if !ool
      } //if a1-a2 and a3-a4 bonds exist
    } //if nadd==2, loop over a1,a2,a3,a4
  } //loop i over nbf
#endif


 //handle zero breaks
  if (nbf==0 && wsh==-1)
  {

    if (nadd==1)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    {
      if (!ic1.bond_exists(a1,a2))
      {
        int ool = 0;
        for (int j=0;j<natoms1;j++) coordn[j] = ic1.coordn[j];
        for (int j=0;j<natoms1;j++) moving[j] = 0;

        coordn[a1]++;
        coordn[a2]++;
        adds[0] = a1;
        adds[1] = a2;
        moving[a1] = 1;
        moving[a2] = 1;
        //printf(" a1: %2i a2: %2i \n",a1+1,a2+1);

        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn);
        if (!ool && nr>1)
          ool = two_frags(a1,a2,ic1);

        int tm_in = 0;
        if (ic1.isTM(a1)) tm_in++; if (ic1.isTM(a2)) tm_in++;
        if (!ool && tm_in)
        {
          nfound += tm_angle_drive(nr,nadd,adds,0,NULL,ic1tm,nangles,angles,anglev,ntor,tor,torv,active1);
        }
        else if (!ool)
        {
          int ga_pass = ml_eval(ic1,nadd,adds,0,NULL);
          if (ga_pass)
          {
            if (nr>1)
            {
              align1.add_align(nadd,adds);
              write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
            }
            else
              write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
            write_ISOMER(niso++,nadd,adds,0,NULL);
            nfound++;
          }
        } //if !ool
      } //if a1-a2 bond exists
    } //if nadd==1, loop over a1,a2


   //CPMZ check me
    if (nadd==2)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    for (int a3=0;a3<=a1;a3++)
    if (active1[a3]==1)
    for (int a4=0;a4<a3;a4++)
    if (active1[a4]==1)
    if (!(a1==a3 && a4>=a2))
    {
      if (!ic1.bond_exists(a1,a2) && !ic1.bond_exists(a3,a4))
      {
        int ool = 0;
        for (int j=0;j<natoms1;j++) coordn[j] = ic1.coordn[j];
        for (int j=0;j<natoms1;j++) moving[j] = 0;

        coordn[a1]++;
        coordn[a2]++;
        coordn[a3]++;
        coordn[a4]++;
        adds[0] = a1;
        adds[1] = a2;
        adds[2] = a3;
        adds[3] = a4;
        moving[a1] = 1;
        moving[a2] = 1;
        moving[a3] = 1;
        moving[a4] = 1;
        //printf(" a1: %2i a2: %2i a3: %2i a4 %2i \n",a1+1,a2+1,a3+1,a4+1);

        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn) + hh_bond(a3,a4,ic1.anumbers,coordn);
        if (!ool && nr>1)
          ool = two_frags(a1,a2,a3,a4,ic1);

        int tm_in = 0;
        if (ic1.isTM(a1)) tm_in++; if (ic1.isTM(a2)) tm_in++;
        if (ic1.isTM(a3)) tm_in++; if (ic1.isTM(a4)) tm_in++;
        if (!ool && tm_in)
        {
          nfound += tm_angle_drive(nr,nadd,adds,0,NULL,ic1tm,nangles,angles,anglev,ntor,tor,torv,active1);
        }
        else if (!ool)
        {
          int ga_pass = ml_eval(ic1,nadd,adds,0,NULL);
          if (ga_pass)
          {
            if (nr>1)
            {
              align1.add_align(nadd,adds);
              write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
            }
            else
              write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
            write_ISOMER(niso++,nadd,adds,0,NULL);
            nfound++;
          }
        }
      }
    } //if nadd==2, loop over a1,a2,a3,a4

  } //if nbrks==0

  delete [] satoms;
  delete [] brksp;
  delete [] adds;
  delete [] coordn;
  delete [] coordnb;
  delete [] moving;
  delete [] movingb;

  delete [] nangles;
  delete [] ntor;
  for (int i=0;i<30;i++)
    delete [] angles[i];
  for (int i=0;i<30;i++)
    delete [] tor[i];
  delete [] angles;
  delete [] tor;


  printf("   found %2i isomers, %2i skipped due to R-H-H bonds \n",nfound,nhh);
  return nfound;
} 


int ZStruct::get_combo(int nadd, int nbf, int nbrks, int* brks, ICoord ic1, int nr, int* active1, int wsh)
{
  if (nadd+nbf<1) return 0;

  printf("  in get_combo: %i/%i \n",nadd,nbrks);

  int maxchg = MAX_CHG;
  int minchg = MIN_CHG;
  if (nadd+nbrks > maxchg)
  {
    printf("  connectivity change limit exceeded (%i/%i) \n",nadd+nbrks,maxchg);
    return 0;
  }
  if (nadd+nbrks<minchg)
  {
    printf("   connectivity change too small (%i/%i) \n",nadd+nbrks,minchg);
    return 0;
  }

  //printf(" getting combo with %i atoms \n",ic1.natoms); fflush(stdout);
  //ic1.print_ic();

  int natoms1 = ic1.natoms;
  int* coordn = new int[natoms1];
  int* coordnb = new int[natoms1]; //after breaks are set

  int nfound = 0;
  int* adds = new int[nadd*2+1];

  int* moving = new int[natoms1];
  int* movingb = new int[natoms1];
  nhh = 0;

  int* satoms = new int[3];


#if 1
  for (int i=0;i<nbf;i++)
  {
    for (int j=0;j<natoms1;j++) coordnb[j] = ic1.coordn[j];
    for (int j=0;j<natoms1;j++) movingb[j] = 0;

    for (int j=0;j<nbrks;j++)
    {
      int b1 = brks[2*nbrks*i+2*j+0];
      int b2 = brks[2*nbrks*i+2*j+1];
      coordnb[b1]--;
      coordnb[b2]--;
      movingb[b1] = 1;
      movingb[b2] = 1;
      //printf(" b1: %2i b2: %2i \n",b1+1,b2+1);
    }

    if (nadd==0 && nr==1 && wsh==-1)
    {
      int ool = 0;
      for (int j=0;j<natoms1;j++)
      if (active1[j]!=1 && movingb[j])
        ool = 1;
      if (!ool)
      for (int j=0;j<natoms1;j++)
      if (active1[j]==1 && movingb[j])
      {
        //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
        if (coordnb[j]<climit_l[ic1.anumbers[j]] || coordnb[j]>climit_h[ic1.anumbers[j]])
          ool = 1;
      }
      if (!ool)
      {
        int ga_pass = ml_eval(ic1,0,NULL,nbrks,&brks[2*nbrks*i]);
        if (ga_pass)
        {
          if (nr>1) //should never be here
          {
            align1.align_zero();
            write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
          }
          else
            write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
          write_ISOMER(niso++,0,NULL,nbrks,&brks[2*nbrks*i]);
          nfound++;
        }
      } //if !ool
    } //if nadd==0 && nr==1


    if (nadd==1)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    {
      int ool = 0;
      for (int j=0;j<natoms1;j++) coordn[j] = coordnb[j];
      for (int j=0;j<natoms1;j++) moving[j] = movingb[j];

      if (!ic1.bond_exists(a1,a2))
      {
        coordn[a1]++;
        coordn[a2]++;
        adds[0] = a1;
        adds[1] = a2;
        moving[a1] = 1;
        moving[a2] = 1;
        //printf(" a1: %2i a2: %2i \n",a1+1,a2+1);

        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn);
        if (!ool && nr>1)
          ool = two_frags(a1,a2,ic1);

        if (!ool)
        {
          if (wsh==-1)
          {
            int ga_pass = ml_eval(ic1,nadd,adds,nbrks,&brks[2*nbrks*i]);
            if (ga_pass)
            {
              if (nr>1)
              {
                align1.add_align(nadd,adds);
                write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
              }
              else
                write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
              write_ISOMER(niso++,nadd,adds,nbrks,&brks[2*nbrks*i]);
              nfound++;
            }
          }
          else if (wsh>-1)
          {
            if (nr>1)
              align1.add_align(nadd,adds);
            int wsh1 = is_shuttle(ic1,nadd,adds,nbrks,&brks[2*nbrks*i],satoms);
            if (wsh==wsh1)
            {
              nfound += get_shuttle(ic1,nr,wsh,nadd,adds,nbrks,&brks[2*nbrks*i],satoms);
            }
          }

        } //if !ool
      } //if a1-a2 bond exists
    } //if nadd1==1, loop over a1,a2


   //CPMZ check me
    if (nadd==2)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    for (int a3=0;a3<=a1;a3++)
    if (active1[a3]==1)
    for (int a4=0;a4<a3;a4++)
    if (active1[a4]==1)
    if (!(a1==a3 && a4>=a2))
    {
      int ool = 0;
      for (int j=0;j<natoms1;j++) coordn[j] = coordnb[j];
      for (int j=0;j<natoms1;j++) moving[j] = movingb[j];

      if (!ic1.bond_exists(a1,a2) && !ic1.bond_exists(a3,a4))
      {
        coordn[a1]++;
        coordn[a2]++;
        coordn[a3]++;
        coordn[a4]++;
        adds[0] = a1;
        adds[1] = a2;
        adds[2] = a3;
        adds[3] = a4;
        moving[a1] = 1;
        moving[a2] = 1;
        moving[a3] = 1;
        moving[a4] = 1;
        //printf(" a1: %2i a2: %2i a3: %2i a4 %2i \n",a1+1,a2+1,a3+1,a4+1);

       //not H-transfer version, see above
        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn) + hh_bond(a3,a4,ic1.anumbers,coordn);
        if (!ool && nr>1)
          ool = two_frags(a1,a2,a3,a4,ic1);

        if (!ool)
        {
          if (wsh==-1)
          {
            int ga_pass = ml_eval(ic1,nadd,adds,nbrks,&brks[2*nbrks*i]);
            if (ga_pass)
            {
              if (nr>1)
              {
                align1.add_align(nadd,adds);
                write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
              }
              else
                write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
              write_ISOMER(niso++,nadd,adds,nbrks,&brks[2*nbrks*i]);
              nfound++;
            }
          } //if !shuttle
          else if (wsh>-1)
          {
            if (nr>1)
              align1.add_align(nadd,adds);
            int wsh1 = is_shuttle(ic1,nadd,adds,nbrks,&brks[2*nbrks*i],satoms);
            if (wsh==wsh1)
            {
              nfound += get_shuttle(ic1,nr,wsh,nadd,adds,nbrks,&brks[2*nbrks*i],satoms);
            }
          }

        } //if !ool
      } //if a1-a2 and a3-a4 bonds exist
    } //if nadd==2, loop over a1,a2,a3,a4
  } //loop i over nbf
#endif


 //handle zero breaks
  if (nbf==0 && wsh==-1)
  {

    if (nadd==1)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    {
      if (!ic1.bond_exists(a1,a2))
      {
        int ool = 0;
        for (int j=0;j<natoms1;j++) coordn[j] = ic1.coordn[j];
        for (int j=0;j<natoms1;j++) moving[j] = 0;

        coordn[a1]++;
        coordn[a2]++;
        adds[0] = a1;
        adds[1] = a2;
        moving[a1] = 1;
        moving[a2] = 1;
        //printf(" a1: %2i a2: %2i \n",a1+1,a2+1);

        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn);
        if (!ool && nr>1)
          ool = two_frags(a1,a2,ic1);

        if (!ool)
        {
          int ga_pass = ml_eval(ic1,nadd,adds,0,NULL);
          if (ga_pass)
          {
            if (nr>1)
            {
              align1.add_align(nadd,adds);
              write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
            }
            else
              write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
            write_ISOMER(niso++,nadd,adds,0,NULL);
            nfound++;
          }
        } //if !ool
      } //if a1-a2 bond exists
    } //if nadd==1, loop over a1,a2


   //CPMZ check me
    if (nadd==2)
    for (int a1=0;a1<natoms1;a1++)
    if (active1[a1]==1)
    for (int a2=0;a2<a1;a2++)
    if (active1[a2]==1)
    for (int a3=0;a3<=a1;a3++)
    if (active1[a3]==1)
    for (int a4=0;a4<a3;a4++)
    if (active1[a4]==1)
    if (!(a1==a3 && a4>=a2))
    {
      if (!ic1.bond_exists(a1,a2) && !ic1.bond_exists(a3,a4))
      {
        int ool = 0;
        for (int j=0;j<natoms1;j++) coordn[j] = ic1.coordn[j];
        for (int j=0;j<natoms1;j++) moving[j] = 0;

        coordn[a1]++;
        coordn[a2]++;
        coordn[a3]++;
        coordn[a4]++;
        adds[0] = a1;
        adds[1] = a2;
        adds[2] = a3;
        adds[3] = a4;
        moving[a1] = 1;
        moving[a2] = 1;
        moving[a3] = 1;
        moving[a4] = 1;
        //printf(" a1: %2i a2: %2i a3: %2i a4 %2i \n",a1+1,a2+1,a3+1,a4+1);

        for (int j=0;j<natoms1;j++)
        if (active1[j]!=1 && moving[j])
          ool = 1;
        if (!ool)
        for (int j=0;j<natoms1;j++)
        if (active1[j]==1 && moving[j])
        {
          //printf(" coordn[j]: %i climit_l: %i climit_h: %i \n",coordn[j],climit_l[ic1.anumbers[j]],climit_h[ic1.anumbers[j]]);
          if (coordn[j]<climit_l[ic1.anumbers[j]] || coordn[j]>climit_h[ic1.anumbers[j]])
            ool = 1;
        }
        if (!ool)
          ool = hh_bond(a1,a2,ic1.anumbers,coordn) + hh_bond(a3,a4,ic1.anumbers,coordn);
        if (!ool && nr>1)
          ool = two_frags(a1,a2,a3,a4,ic1);

        if (!ool)
        {
          int ga_pass = ml_eval(ic1,nadd,adds,0,NULL);
          if (ga_pass)
          {
            if (nr>1)
            {
              align1.add_align(nadd,adds);
              write_initial_xyz(niso,natoms1,ic1.anames,align1.xyza,ic1.q1);
            }
            else
              write_initial_xyz(niso,natoms1,ic1.anames,ic1.coords,ic1.q1);
            write_ISOMER(niso++,nadd,adds,0,NULL);
            nfound++;
          }
        }
      }
    } //if nadd==2, loop over a1,a2,a3,a4

  } //if nbrks==0

  delete [] satoms;
  delete [] adds;
  delete [] coordn;
  delete [] coordnb;
  delete [] moving;
  delete [] movingb;

  printf("  found %2i isomers, %2i skipped due to R-H-H bonds \n",nfound,nhh);
  return nfound;
} 



int ZStruct::get_breaks(int nbrks, int* brks, ICoord ic1, int* active1)
{
  //printf("\n in get_breaks: %i \n",nbrks);

  if (nbrks<1) return 0;

 //# of possible break combinations
  int nbt = ic1.nbonds;
  for (int i=1;i<nbrks;i++)
    nbt *= (ic1.nbonds-i);
  nbt = nbt/fact(nbrks);

  for (int i=0;i<2*nbt*nbrks;i++)
    brks[i] = -1;


  int nfound = 0;

  if (nbrks==1)
  for (int i=0;i<ic1.nbonds;i++)
  {
    int a1 = ic1.bonds[i][0];
    int a2 = ic1.bonds[i][1];
    if (active1[a1]==1 && active1[a2]==1)
    {
      brks[2*nfound+0] = a1;
      brks[2*nfound+1] = a2;
      nfound++;
    }
  }

  if (nbrks==2)
  for (int i=0;i<ic1.nbonds;i++)
  for (int j=0;j<i;j++)
  {
    int a1 = ic1.bonds[i][0];
    int a2 = ic1.bonds[i][1];
    int a3 = ic1.bonds[j][0];
    int a4 = ic1.bonds[j][1];
    if (active1[a1]==1 && active1[a2]==1
     && active1[a3]==1 && active1[a4]==1)
    {
      brks[4*nfound+0] = a1;
      brks[4*nfound+1] = a2;
      brks[4*nfound+2] = a3;
      brks[4*nfound+3] = a4;
      nfound++;
    }
  }

  if (nbrks==3)
  for (int i=0;i<ic1.nbonds;i++)
  for (int j=0;j<i;j++)
  for (int k=0;k<j;k++)
  {
    int a1 = ic1.bonds[i][0];
    int a2 = ic1.bonds[i][1];
    int a3 = ic1.bonds[j][0];
    int a4 = ic1.bonds[j][1];
    int a5 = ic1.bonds[k][0];
    int a6 = ic1.bonds[k][1];
    if (active1[a1]==1 && active1[a2]==1 
     && active1[a3]==1 && active1[a4]==1
     && active1[a5]==1 && active1[a6]==1)
    {
      brks[6*nfound+0] = a1;
      brks[6*nfound+1] = a2;
      brks[6*nfound+2] = a3;
      brks[6*nfound+3] = a4;
      brks[6*nfound+4] = a5;
      brks[6*nfound+5] = a6;
      nfound++;
    }
  }

  if (nbrks==4)
  for (int i=0;i<ic1.nbonds;i++)
  for (int j=0;j<i;j++)
  for (int k=0;k<j;k++)
  for (int l=0;l<k;l++)
  {
    int a1 = ic1.bonds[i][0];
    int a2 = ic1.bonds[i][1];
    int a3 = ic1.bonds[j][0];
    int a4 = ic1.bonds[j][1];
    int a5 = ic1.bonds[k][0];
    int a6 = ic1.bonds[k][1];
    int a7 = ic1.bonds[l][0];
    int a8 = ic1.bonds[l][1];
    if (active1[a1]==1 && active1[a2]==1
     && active1[a3]==1 && active1[a4]==1
     && active1[a5]==1 && active1[a6]==1
     && active1[a7]==1 && active1[a8]==1)
    {
      brks[8*nfound+0] = a1;
      brks[8*nfound+1] = a2;
      brks[8*nfound+2] = a3;
      brks[8*nfound+3] = a4;
      brks[8*nfound+4] = a5;
      brks[8*nfound+5] = a6;
      brks[8*nfound+6] = a7;
      brks[8*nfound+7] = a8;
      nfound++;
    }
  }

#if 1
  printf("  breaks found (%i):",nfound);
  for (int i=0;i<nfound;i++)
  {
    for (int j=0;j<nbrks;j++)
      printf(" %i-%i",brks[2*nbrks*i+2*j+0]+1,brks[2*nbrks*i+2*j+1]+1);
  }
  printf("\n");
#endif

  return nfound;
}


int ZStruct::get_shuttle(ICoord ic1, int nr, int wsh, int nadd1, int* add1, int nbrks1, int* brks1, int* s1)
{
  int found = 0;

  //printf("  in get_shuttle \n");
  //printf("   atoms shuttle: H-O: %i-%i destination: %i \n",s1[0]+1,s1[1]+1,s1[2]+1);

  int natoms1 = ic1.natoms;
  int natomsn = icshuttle[wsh].natoms;
  int natoms2 = natoms1 + natomsn;
  int nadd2 = nadd1 + 1;
  int nbrks2 = nbrks1 + 1;
  int* add2 = new int[2*nadd2];
  int* brks2 = new int[2*nbrks2];
  int* add3 = new int[4];
  int nadd3 = 2;
  double* xyz2 = new double[3*natoms2];
  string* anames2 = new string[natoms2];
  int* anumbers2 = new int[natoms2];
 //CPMZ fix me
  int q1 = 0;

  int wshr = wshuttle[wsh];
  double* qrc2 = new double[natoms2];
  for (int i=0;i<natoms1;i++) 
    qrc2[i] = ic1.q[i];
  for (int i=0;i<natomsn;i++)
    qrc2[natoms1+i] = qr[wshr][i];

  //printf("   natoms2: %i nadd2: %i nbrks2: %i \n",natoms2,nadd2,nbrks2);

#if 0
  for (int i=0;i<nadd1;i++)
    printf("   adding: %i %i \n",add1[2*i+0]+1,add1[2*i+1]+1);
  for (int i=0;i<nbrks1;i++)
    printf("   breaking: %i %i \n",brks1[2*i+0]+1,brks1[2*i+1]+1);
#endif

  //fill add2 and brks2 vectors
  int naddf = 0;
  int a1 = s1[0];
  int a2 = natoms1 + shuttles[wsh][2];
  int a3 = s1[2];
  int a4 = natoms1 + shuttles[wsh][0];
  for (int i=0;i<nadd1;i++)
  {
    if (!(add1[2*i+0]==a1 && add1[2*i+1]==a3)
        && !(add1[2*i+0]==a3 && add1[2*i+1]==a1))
    {
      add2[2*naddf+0] = add1[2*i+0];
      add2[2*naddf+1] = add1[2*i+1];
      naddf++;
    }
  }
  //printf("   naddf: %i nadd1: %i nadd2: %i \n",naddf,nadd1,nadd2);
  add2[2*naddf+0] = a1;
  add2[2*naddf+1] = a2;
  add2[2*naddf+2] = a3;
  add2[2*naddf+3] = a4;
  for (int i=0;i<2*nbrks1;i++)
    brks2[i] = brks1[i];
  brks2[2*nbrks1+0] = natoms1 + shuttles[wsh][0];
  brks2[2*nbrks1+1] = natoms1 + shuttles[wsh][1];

  add3[0] = a1; add3[1] = a2; add3[2] = a3; add3[3] = a4;

#if 0
  printf("   after shuttle change \n");
  for (int i=0;i<nadd2;i++)
    printf("   adding: %i %i \n",add2[2*i+0]+1,add2[2*i+1]+1);
  for (int i=0;i<nbrks2;i++)
    printf("   breaking: %i %i \n",brks2[2*i+0]+1,brks2[2*i+1]+1);
#endif
//  fflush(stdout);

//  printf(" ic1.natoms: %i align1.natoms1: %i align1.natoms2: %i \n",natoms1,align1.natoms1,align1.natoms2);

  //align the three
  for (int i=0;i<3*natoms1;i++)
//    xyz2[i] = align1.xyza[i]; //uninitialized
    xyz2[i] = ic1.coords[i];
  for (int i=0;i<3*natomsn;i++)
    xyz2[3*natoms1+i] = icshuttle[wsh].coords[i]-20.;
  for (int i=0;i<natoms1;i++)
  {
    anames2[i] = ic1.anames[i];
    anumbers2[i] = ic1.anumbers[i];
  }
  for (int i=0;i<natomsn;i++)
  {
    anames2[natoms1+i] = icshuttle[wsh].anames[i];
    anumbers2[natoms1+i] = icshuttle[wsh].anumbers[i];
  }

//  printf("   before alignment: \n"); 
//  print_xyz_gen(natoms2,anames2,xyz2);


  ICoord ic2;
  ic2.init(natoms2,anames2,anumbers2,xyz2);
  ic2.q = qrc2;

  int ga_pass = ml_eval(ic2,nadd2,add2,nbrks2,brks2);
  if (ga_pass)
  {
    if (nr==2)
    {
      align1.add_third(icshuttle[wsh].natoms,icshuttle[wsh].anames,icshuttle[wsh].anumbers,icshuttle[wsh].coords);
      //align1.add_align already called
      align1.shuttle_align(nadd3,add3);
      write_initial_xyz(niso,natoms2,anames2,align1.xyza3,q1);
      //print_xyz_gen(natoms2,anames2,align1.xyza3);
    }
    else
    {
      align1.init(ic1.natoms,ic1.anames,ic1.anumbers,ic1.coords,icshuttle[wsh].natoms,icshuttle[wsh].anames,icshuttle[wsh].anumbers,icshuttle[wsh].coords);
      align1.add_align(nadd2,add2);
      write_initial_xyz(niso,natoms2,anames2,align1.xyza,q1);
      //print_xyz_gen(natoms2,anames2,align1.xyza);
      align1.freemem();
    }

    write_ISOMER(niso++,nadd2,add2,nbrks2,brks2);
 
    found = 1;
  } //if ga_pass


  ic2.freemem();

  delete [] add2;
  delete [] brks2;
  delete [] add3;
  delete [] xyz2;
  delete [] anames2;
  delete [] anumbers2;
  delete [] qrc2;

  return found;
}


int ZStruct::is_shuttle(ICoord ic1, int nadd1, int* add1, int nbrks1, int* brks1, int* s1)
{
 //currently checks for 3 and 4 membered ring shuttles

  int isShuttle = 0;

  if (nadd1<1 || nbrks1<1 || nshuttle==0)
    return -1;

#if 0
  printf(" in is_shuttle, shuttle type: %i \n",shuttles[0][3]);
  for (int i=0;i<nadd1;i++)
    printf("   adding: %i %i \n",add1[2*i+0]+1,add1[2*i+1]+1);
  for (int i=0;i<nbrks1;i++)
    printf("   breaking: %i %i \n",brks1[2*i+0]+1,brks1[2*i+1]+1);
#endif

#if 1
//takes first shuttle that is found
 //priority to brks order (due to get_h_transfer)
  int wsh = -1;
  for (int n=0;n<nshuttle;n++)
  if (wsh==-1)
  {
    int e1a = shuttles[n][3]; //element shuttled

    for (int i=0;i<nbrks1;i++)
    if (wsh==-1)
    {
      int b1 = brks1[2*i+0];
      int b2 = brks1[2*i+1];
      int e1 = ic1.anumbers[b1];
      int e2 = ic1.anumbers[b2];

      if (e1==e1a || e2==e1a)
      for (int j=0;j<nadd1;j++)
      {
        int a1 = add1[2*j+0];
        int a2 = add1[2*j+1];

        //printf("    e1/e2: %i %i   a1/a2: %i %i b1/b2: %i %i \n",e1,e2,a1,a2,b1,b2);
        if (b1==a1 && e1a==e1)
        {
          //printf("   potential shuttle atom: %i (element: %i) \n",a1+1,e1);
          wsh = n;
          s1[0] = a1;
          s1[1] = b2;
          s1[2] = a2;
          isShuttle++;
          break;
        }
        else if (b1==a2 && e1a==e1)
        {
          //printf("   potential shuttle atom: %i (element: %i) \n",a2+1,e2);
          wsh = n;
          s1[0] = a2;
          s1[1] = b2;
          s1[2] = a1;
          isShuttle++;
          break;
        }
        else if (b2==a1 && e1a==e2)
        {
          //printf("   potential shuttle atom: %i (element: %i) \n",a1+1,e1);
          wsh = n;
          s1[0] = a1;
          s1[1] = b1;
          s1[2] = a2;
          isShuttle++;
          break;
        }
        else if (b2==a2 && e1a==e2)
        {
          //printf("   potential shuttle atom: %i (element: %i) \n",a2+1,e2);
          wsh = n;
          s1[0] = a2;
          s1[1] = b1;
          s1[2] = a1;
          isShuttle++;
          break;
        }
      } //loop j over nbrks1
    } //loop i over nadd1

  } //loop n over nshuttles
#else
//takes first shuttle that is found
 //old code, priority to add list
  int wsh = -1;
  for (int n=0;n<nshuttle;n++)
  if (wsh==-1)
  {
    int e1a = shuttles[n][3]; //element shuttled

    for (int i=0;i<nadd1;i++)
    if (wsh==-1)
    {
      int a1 = add1[2*i+0];
      int a2 = add1[2*i+1];
      int e1 = ic1.anumbers[a1];
      int e2 = ic1.anumbers[a2];

      if (e1==e1a || e2==e1a)
      for (int j=0;j<nbrks1;j++)
      {
        int b1 = brks1[2*j+0];
        int b2 = brks1[2*j+1];

        if (b1==a1 && e1a==e1)
        {
          //printf(" potential shuttle atom: %i (element: %i) \n",a1+1,e1);
          wsh = n;
          s1[0] = a1;
          s1[1] = b2;
          s1[2] = a2;
          isShuttle++;
          break;
        }
        else if (b1==a2 && e1a==e1)
        {
          //printf(" potential shuttle atom: %i (element: %i) \n",a2+1,e2);
          wsh = n;
          s1[0] = a2;
          s1[1] = b2;
          s1[2] = a1;
          isShuttle++;
          break;
        }
        else if (b2==a1 && e1a==e2)
        {
          //printf(" potential shuttle atom: %i (element: %i) \n",a1+1,e1);
          wsh = n;
          s1[0] = a1;
          s1[1] = b1;
          s1[2] = a2;
          isShuttle++;
          break;
        }
        else if (b2==a2 && e1a==e2)
        {
          //printf(" potential shuttle atom: %i (element: %i) \n",a2+1,e2);
          wsh = n;
          s1[0] = a2;
          s1[1] = b1;
          s1[2] = a1;
          isShuttle++;
          break;
        }
      } //loop j over nbrks1
    } //loop i over nadd1

  } //loop n over nshuttles
#endif

 //check if H is involved in > 1 add
  int nhadd = 0;
  if (isShuttle)
  for (int i=0;i<2*nadd1;i++)
  if (add1[i]==s1[0])
    nhadd++; 
  if (nhadd>1)
  {
    //printf("  skipping this shuttle, H involved in add %i times \n",nhadd);
    isShuttle = 0;
  }

 //check for ring
  int isRing = 0;

  // 3 members connected and changing
  if (isShuttle && ic1.bond_exists(s1[1],s1[2]))
  {
    //printf("  ring: %i %i %i \n",s1[0]+1,s1[1]+1,s1[2]+1);
    isRing = 3;
  }

  // 4th member in between endpoints 
  if (!isRing && isShuttle && nbrks1==1 && nadd1==1)
  {
    for (int i=0;i<ic1.natoms;i++)
    if (ic1.bond_exists(s1[1],i))
    if (ic1.bond_exists(i,s1[2]))
    {
      //printf("  possible ring: %i %i %i %i \n",s1[0]+1,s1[1]+1,s1[2]+1,i+1);
      isRing = 4;
      break;
    }
  }

  //
  int r1,r1b,r1a;
  if (nbrks1>1 && isShuttle && !isRing)
  for (int i=0;i<nbrks1;i++)
  {
    int b1 = brks1[2*i+0];
    int b2 = brks1[2*i+1];
    r1 = -1;
    if (b1==s1[2])
    {
      r1 = b2;
      r1b = b1;
    }
    if (b2==s1[2])
    {
      r1 = b1;
      r1b = b2;
    }
    if (r1>-1)
    {
      //printf("  possible (brk) ring atoms: %i %i %i %i \n",s1[1]+1,s1[0]+1,s1[2]+1,r1+1);
      if (ic1.bond_exists(r1,s1[1])) 
      {
        isRing = 4;
      }
    }
  } //loop i over nbrks1


  if (nadd1>1 && isShuttle && !isRing)
  for (int i=0;i<nadd1;i++)
  {
    int a1 = add1[2*i+0];
    int a2 = add1[2*i+1];
    r1 = -1;

    if (a1==s1[1])
    {
      r1 = a2;
      r1a = a1;
    }
    if (a2==s1[1])
    {
      r1 = a1;
      r1a = a2;
    }

    if (r1>-1)
    {
      //printf("  possible (add) ring atoms: %i %i %i %i \n",s1[1]+1,s1[0]+1,s1[2]+1,r1+1);
      if (ic1.bond_exists(r1,s1[2]))
      {
        isRing = 4;
      }
    }
  } //loop i over nadd1


//  if (isRing) printf(" %i-membered ring! \n",isRing);
  if (isRing>RINGSHUTTLE || isRing==0)
    isShuttle = 0;

  //if (isShuttle)
  //  printf("  isShuttle: %i wsh: %i ring size: %i \n",isShuttle,wsh,isRing);
  if (isShuttle>1) 
    printf("  WARNING: two H shuttled? \n");


  if (isShuttle==0) wsh = -1;
  return wsh;
}


int ZStruct::two_frags(int a1, int a2, ICoord ic1)
{
 //returns 1 if fragments are equal
  int result = 1;

  if (ic1.frags[a1]!=ic1.frags[a2])
    result = 0;

  //printf(" in two_frags: %i %i result: %i \n",a1,a2,result);

  return result;
}

int ZStruct::two_frags(int a1, int a2, int a3, int a4, ICoord ic1)
{
 //returns 1 if fragments are equal
  int result = 1;

  if (ic1.frags[a1]!=ic1.frags[a2])
    result = 0;
  else if (ic1.frags[a3]!=ic1.frags[a4])
    result = 0;

  //printf(" in two_frags: %i %i %i %i result: %i \n",a1,a2,a3,a4,result);

  return result;
}


int ZStruct::pair_allowed(int i, int j, ICoord* icr)
{
  int allowed = 1;

  int natoms1 = icr[i].natoms;
  int natoms2 = icr[j].natoms;

  if (natoms1+natoms2>atomlimit)
    allowed = 0;

  //printf(" pair_allowed: %i %i total: %i allowed? %i \n",natoms1,natoms2,natoms1+natoms2,allowed);

  return allowed;
}

int ZStruct::hh_bond(int a1, int a2, int* anumbers, int* coordn)
{
  int ool = 0;

 //2 coord H-H bond
  if ((anumbers[a1]==1 && anumbers[a2]==1) && (coordn[a1]>1 || coordn[a2]>1))
    ool = 1;

//  printf("     hhbond: %i \n",ool);
  if (ool)
    nhh++;

  return ool;
}


void ZStruct::init()
{
  struct stat sts;
  if (stat("scratch",&sts) == -1)
  {
    printf("\n ERROR: no scratch/ directory present, exiting \n");
    exit(1);
  }

  ga_avail = 0;
  temperature = 298.15;
  pthresh = 0.5;
  atomlimit = 1000;

  climit_l = NULL;
  climit_h = NULL;
  activea = 0;
  active = NULL;

  nreact = 0;
  nreacta = 0;
  reactnum = NULL;
  reactrun = NULL;
  npair = 0;
  pairs = NULL;
  pairsnum = NULL;
  wpaira = 0;
  rsaved = NULL;
  wpair = NULL;
  wshpair = NULL;
  upair = NULL;

  natoms = NULL;
  anames = NULL;
  anumbers = NULL;
  xyzr = NULL;
  xyzp = NULL;
  xyzt = NULL;
  qrc = NULL;

  pair_react = NULL;
  natomsr = NULL;
  qr = NULL;
  qp = NULL; //not same size as qr

  nshuttle = 0;
  wshuttle = NULL;
  shuttles = NULL;
  icshuttle = NULL;

  elista = 0;
  elistr = NULL;
  elistrp = NULL;
  elistc = NULL;
  elistcmin = NULL;
  elistcminsh = NULL;
  elistts = NULL;
  elistp = NULL;
  elistref = NULL;
  Ea = NULL;
  Erxn = NULL;

  nfragb = 0;
  fragb = NULL;

  nbo_reacts = NULL;
  rxns1.init(1); //level one rxn analysis
  rxns2.init(1); 

  align1.inited = 0;

  niso = 0;
  nisoa = 0;

  nhh = 0;

  return;
}

void ZStruct::dealloc()
{
  if (climit_l!=NULL) delete [] climit_l;
  if (climit_h!=NULL) delete [] climit_h;

  if (elistr!=NULL)  delete [] elistr;
  if (elistrp!=NULL) delete [] elistrp;
  if (elistc!=NULL)  delete [] elistc;
  if (elistts!=NULL) delete [] elistts;
  if (elistp!=NULL)  delete [] elistp;

  return;
}


void ZStruct::fill_elistrp()
{
  //printf(" in fill_elistrp, npair: %i \n",npair);
  if (elistrp!=NULL) delete [] elistrp;
  elistrp = new double[npair];

  int nfound = 0;
#if 0
  for (int i=0;i<nreact;i++)
  for (int j=0;j<=i;j++)
    elistrp[nfound++] = elistr[i] + elistr[j];
#endif
  for (int i=0;i<npair;i++)
  {
    int wr1 = pairs[2*i+0];
    int wr2 = pairs[2*i+1];
    elistrp[i] = elistr[wr1] + elistr[wr2];
  }

  if (nfound>npair) printf(" WARNING: nfound>npair in fill_elistrp \n");

#if 0
  for (int i=0;i<npair;i++)
    printf(" pair %2i total E: %8.6f \n",i,elistrp[i]);
#endif

  return;
}

void ZStruct::calc_Ea_Erxn()
{
  printf(" in calc_Ea_Erxn() \n"); fflush(stdout);

  if (Ea!=NULL) delete [] Ea;
  if (Erxn!=NULL) delete [] Erxn;

  Ea = new double[niso];
  Erxn = new double[niso];
  for (int i=0;i<niso;i++) Ea[i] = 111.1;
  for (int i=0;i<niso;i++) Erxn[i] = 0.;

  for (int i=0;i<niso;i++)
  if (rsaved[i])
  {
    //printf(" iso: %3i \n",i); fflush(stdout);
    int wp = wpair[i];
    double ets; 
    double esh = 0.;
    int wsh = wshpair[i];
    if (wsh>-1) esh = elistr[wshuttle[wsh]];
    //printf(" esh: %8.6f ",esh); //printf(" wp: %i ",wp);

    if (wp>=0 && wsh==-1)
      ets = (elistts[i] - elistcmin[wp])*627.5;
//      ets = (elistts[i] - elistrp[wp] - esh)*627.5;
    else if (wp>=0 && wsh>-1)
      ets = (elistts[i] - elistcminsh[wp])*627.5;
    else
      ets = (elistts[i] - elistref[i])*627.5; //was elistc

   //CPMZ fix me.. don't use c (?)
    double erxn;
    if (wp>=0) 
      erxn = (elistp[i] - elistrp[wp] - esh)*627.5;
    else 
      erxn = (elistp[i] - elistref[i])*627.5; //was elistc

    Ea[i] = ets;
    Erxn[i] = erxn;
  } //loop i over niso

  printf(" done with calc_Ea_Erxn() \n"); fflush(stdout);

  return;
}

