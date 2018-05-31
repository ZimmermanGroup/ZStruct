// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "nbo.h"
#include "constants.h"
#include "zstruct.h"
using namespace std;

//TO DO:
// 1. read in/save molden GTO section (done)
// 2. read in/save hybridization data
// 3. MOPAC data (done-ish)


// NOTE:                          
// bmo_occ doesn't match bmo_ quantities if 3-center is	added


double NBO::get_pol(int a1, int a2, double& ev1)
{
  double pol1 = 0.;
  for (int j=0;j<bmo;j++) //finds bond including atoms a1 and a2
  if ( (bmo_atoms[2*j+0]==a1 && bmo_atoms[2*j+1]==a2)
    || (bmo_atoms[2*j+0]==a2 && bmo_atoms[2*j+1]==a1) )
  {
    if (mo_occ[j]>ev1) //grab highest energy localized orbital
    {
      ev1 = mo_occ[j];
      pol1 = bmo_polar[j];
      if (bmo_atoms[2*j+0]==a1)
        pol1 *= -1.;
    }
  }
  return pol1;
}


void NBO::print_molden_orbs(int norbs, int* olist, string filename)
{
  printf("  printing molden output to %s \n",filename.c_str());

  ofstream molfile;
  molfile.open(filename.c_str());
  molfile.setf(ios::scientific);
  molfile.setf(ios::left);
  molfile << setprecision(8);

  double* diag = mo_occ; //for mopac, this is the "eigenvalue"

  for (int i=0;i<mheadersize;i++)
    molfile << mheader[i] << endl;


  for (int i=0;i<norbs;i++)
  if (olist[i])
  {
    int index = i;

    molfile << "Sym=X" << endl;
    molfile << "Ene= " << diag[index] << endl;
    molfile << "Spin=Alpha" << endl;
    molfile << "Occup=1" << endl;

    for (int j=0;j<nao;j++)
      molfile << "   " << j+1 << "   " << MO[index*nao+j] << endl;
  }
  for (int i=0;i<nao;i++)
  if (!olist[i])
  {
    int index = i;

    molfile << "Sym=X" << endl;
    molfile << "Ene= " << diag[index] << endl;
    molfile << "Spin=Alpha" << endl;
    molfile << "Occup=0" << endl;

    for (int j=0;j<nao;j++)
      molfile << "   " << j+1 << "   " << MO[index*nao+j] << endl;
  }

  molfile << "[6D]" << endl;
  molfile.close();

  return;
}


int NBO::compare_nbo(NBO nbo1, string mfilename, int quiet)
{
  nb = 0;
  if (!hasNBO || !nbo1.hasNBO || nao<1)
  {
    printf("  missing data in compare_nbo \n");
    return -1;
  }

 // printf(" bmo: %2i bmo': %2i \n",bmo,nbo1.bmo);
 // printf(" vmo: %2i vmo': %2i \n",vmo,nbo1.vmo);

  if (blist!=NULL) delete [] blist;
  blist = new int[bmo+1];
  for (int i=0;i<bmo;i++) blist[i] = 0;

  int* foundbmo = new int[bmo+1];
  int* used = new int[nbo1.bmo+1];
  for (int i=0;i<bmo;i++) foundbmo[i] = 0;
  for (int i=0;i<nbo1.bmo;i++) used[i] = 0;
  for (int i=0;i<bmo;i++)
  {
    int found = 0;

    double pol1 = bmo_polar[i];
    int b1 = bmo_atoms[2*i+0];
    int b2 = bmo_atoms[2*i+1];
    for (int j=0;j<nbo1.bmo;j++)
    if (!used[j])
    {
      double pol2 = nbo1.bmo_polar[j];
      int c1 = nbo1.bmo_atoms[2*j+0];
      int c2 = nbo1.bmo_atoms[2*j+1];
      if (c1==b1 && c2==b2) found = 1;
      if (c1==b2 && c2==b1) found = 1;
      if (!found && b2==0) //need to check this block
      {
        if (b1==c1 && pol2<-0.9)
          found = 1;
        else if (b1==c2 && pol2>0.9)
          found = 1;
      }
      if (!found && pol1<-0.9) //need to check this block
      {
        if (b1==c1 && pol2<-0.9)
          found = 1;
        else if (b1==c2 && pol2>0.9)
          found = 1;
      }
      if (found) { used[j] = 1; break; }
    }
    foundbmo[i] = found;
  }

  if (!quiet)
  for (int i=0;i<bmo;i++)
  if (!foundbmo[i])
    printf("   orbital %2i: %2i-%2i \n",i+1,bmo_atoms[2*i+0],bmo_atoms[2*i+1]);

  int* saveorb = new int[nao];
  for (int i=0;i<bmo;i++) blist[i] = 0;
  for (int i=0;i<nao;i++) saveorb[i] = 0;
  for (int i=0;i<bmo;i++)
  if (!foundbmo[i])
  {
    if (bmo_num[i]-1<0 || bmo_num[i]>nao)
    { printf(" ERROR bmo_num out of bounds: %i \n",bmo_num[i]); exit(1); }
    saveorb[bmo_num[i]-1] = 1;
    blist[nb++] = i;
  }
  if (mfilename!="none")
    print_molden_orbs(nao,saveorb,mfilename);

  delete [] foundbmo;
  delete [] saveorb;

  return 0;
}

void NBO::print_nbo()
{
  printf("\n Printing NBO data! \n");
  for (int i=0;i<natoms;i++)
    printf("  q[%2i]: %8.5f \n",i+1,q[i]);
  printf("  nalpha: %2i nbeta: %2i nelec: %2i nao: %2i \n",nalpha,nbeta,nelec,nao);
  printf("\n");
  printf(" now printing NBO bonding analysis \n");

  for (int i=0;i<bmo;i++)
  {
    printf("  BMO[%2i] %2i-%2i occ: %5.3f polarization: %6.3f \n",i+1,bmo_atoms[2*i+0],bmo_atoms[2*i+1],bmo_occ[i],bmo_polar[i]);
  }
  for (int i=0;i<vmo;i++)
  {
    printf("  VMO[%2i] %2i-%2i occ: %5.3f \n",i+1,vmo_atoms[2*i+0],vmo_atoms[2*i+1],vmo_occ[i]);
  }
  printf("\n");
}

void NBO::alloc_nbo()
{
  //printf("  in alloc_nbo for nao: %2i \n",nao);
  int resetnao = 0;
  if (nao<1)
  {
    //printf("  nao reset from: %i \n",nao);
    nao = 10;
    resetnao = 1;
  }
  nao+=1000;
  mo_occ = new double[nao];
  bmo_occ = new double[nao];
  bmo_atoms = new int[2*nao];
  bmo_polar = new double[nao];
  bmo_num = new int[nao];
  vmo_occ = new double[nao];
  vmo_atoms = new int[2*nao];
  vmo_num = new int[nao];
  for (int i=0;i<nao;i++)   bmo_occ[i] = 0;
  for (int i=0;i<2*nao;i++) bmo_atoms[i] = 0;
  for (int i=0;i<nao;i++) bmo_num[i] = 0;
  for (int i=0;i<nao;i++)   vmo_occ[i] = 0;
  for (int i=0;i<2*nao;i++) vmo_atoms[i] = 0;
  for (int i=0;i<nao;i++) vmo_num[i] = 0;
  nao-=1000;

  if (resetnao) nao = 0;

  return;
}

double NBO::distance(int i, int j)
{
  return sqrt((xyz[3*i+0]-xyz[3*j+0])*(xyz[3*i+0]-xyz[3*j+0])+
              (xyz[3*i+1]-xyz[3*j+1])*(xyz[3*i+1]-xyz[3*j+1])+
   	      (xyz[3*i+2]-xyz[3*j+2])*(xyz[3*i+2]-xyz[3*j+2]));
}

void NBO::clean_long_bonds()
{
  //printf("  removing long distance \"bonds\" \n");
  for (int i=0;i<bmo;i++)
  {
    int a1 = bmo_atoms[2*i+0]-1;
    int a2 = bmo_atoms[2*i+1]-1;
    double d1 = 0.; 
    //printf("    a1,a2: %2i %2i anumbers: %i %i \n",a1,a2,anumbers[a1],anumbers[a2]);

    double d0 = 10.;
    if (a1>-1 && a2>-1)
    {
      d0 = getRa(anumbers,a1) + getRa(anumbers,a2);
      d1 = distance(a1,a2);
    }
    //printf("    d0: %8.5f d1: %8.5f \n",d0,d1);
    if (d1 > d0)
    {
      //printf("  long: %2i-%2i dist: %5.2f \n",a1,a2,d1);
      for (int j=i;j<bmo-1;j++)
      {
        int jp = j+1;
        bmo_num[j] = bmo_num[jp];
        bmo_atoms[2*j+0] = bmo_atoms[2*jp+0];
        bmo_atoms[2*j+1] = bmo_atoms[2*jp+1];
        bmo_occ[j] = bmo_occ[jp];
        bmo_polar[j] = bmo_polar[jp];
      }
      bmo--;
    }
  }

  return;
}

int NBO::read_nbo_file(string filename) 
{
  hasNBO = 0;
  bmo = vmo = 0;

#if USE_MOPAC
  int done1 = read_mopac_mo(filename);
  if (done1) hasNBO = 1;
  clean_long_bonds();
  return done1;
  printf(" SHOULDN'T BE HERE \n");
  exit(1);
#endif

  int done = read_mo(filename);
  if (!done) 
  {
    printf(" couldn't read MO's \n");
    return 0;
  }

  alloc_nbo();

  string oname = filename;
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening NBO file: %s \n",oname.c_str()); return 0; }

  string line;
  vector<string> tok_line;
  int cmo = 0;
  while (!output.eof()) 
  { 
    getline(output,line);
   // cout << " RR " << line << endl;

    if (line.find(" beta electrons")!=string::npos)
    {
     // cout << " RR1 " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      nalpha = atoi(tok_line[2].c_str());
      nbeta = atoi(tok_line[5].c_str());
      nelec = nalpha+nbeta;
    }
    if (line.find("Summary of Natural Population Analysis")!=string::npos)
    {
     // cout << " RR1 " << line << endl;
      for (int i=0;i<5;i++)
        getline(output,line);
      for (int i=0;i<natoms;i++)
      {
        getline(output,line);
       // cout << "  RR2 " << line << endl;
        tok_line = StringTools::tokenize(line, " \t");
        q[i] = atof(tok_line[2].c_str());
      }
    }

    if (line.find("Bond orbital/ Coefficients/ Hybrids")!=string::npos)
    {
      while (!output.eof())
      {
       // cout << " RRo " << line << endl;
        getline(output,line);
        if (line.find("BD*")!=string::npos)
        {
          tok_line = StringTools::tokenize(line, "-() \t");
          mo_occ[cmo++] = atof(tok_line[1].c_str());
          vmo_occ[vmo] = atof(tok_line[1].c_str());
          vmo_atoms[2*vmo+0] = atoi(tok_line[5].c_str());
          vmo_atoms[2*vmo+1] = atoi(tok_line[7].c_str());
          vmo_num[vmo] = atoi(tok_line[0].c_str());
          vmo++;
        }
        else if (line.find("RY*")!=string::npos)
        {
          tok_line = StringTools::tokenize(line, "-() \t");
          mo_occ[cmo++] = atof(tok_line[1].c_str());
        }
        else if (line.find("CR")!=string::npos)
        {
          tok_line = StringTools::tokenize(line, "-() \t");
          mo_occ[cmo++] = atof(tok_line[1].c_str());
        }
        else if (line.find("BD")!=string::npos)
        {
          tok_line = StringTools::tokenize(line, "-() \t");
          mo_occ[cmo++] = atof(tok_line[1].c_str());
          bmo_occ[bmo] = atof(tok_line[1].c_str());
          bmo_atoms[2*bmo+0] = atoi(tok_line[5].c_str());
          bmo_atoms[2*bmo+1] = atoi(tok_line[7].c_str());
          bmo_num[bmo] = atoi(tok_line[0].c_str());
          bmo++;
        }
        else if (line.find("LP")!=string::npos)
        {
          tok_line = StringTools::tokenize(line, "-() \t");
          mo_occ[cmo++] = atof(tok_line[1].c_str());
          bmo_occ[bmo] = atof(tok_line[1].c_str());
          bmo_atoms[2*bmo+0] = atoi(tok_line[5].c_str());
          bmo_atoms[2*bmo+1] = 0;
          bmo_num[bmo] = atoi(tok_line[0].c_str());
          bmo++;
        }
        else if (line.find("****************")!=string::npos)
          break;
      }
    }

  }
  output.close();

  hasNBO = 1;

  return 1;
}

void NBO::alloc(int natoms_i)
{ 
  hasNBO = 0;
  bmo = 0;
  vmo = 0;
  natoms = natoms_i;

  anumbers = new int[natoms];
  anames = new string[natoms];
  xyz = new double[3*natoms];

  q = new double[natoms];
  MO = NULL;
  mo_occ = NULL;
  bmo_occ = NULL;
  bmo_atoms = NULL;
  vmo_occ = NULL;
  vmo_atoms = NULL;
  wAO = NULL;
  tAO = NULL;

  return;
}

void NBO::init(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i)
{
  hasNBO = 0;
  bmo = 0;
  vmo = 0;
  nao = 0;
  nalpha = nbeta = nelec = 0;
  natoms = natoms_i;

  anumbers = new int[natoms];
  anames = new string[natoms];
  xyz = new double[3*natoms];

  q = new double[natoms]; 
  spratio = new double[natoms];
  MO = NULL;
  mo_occ = NULL;
  bmo_occ = NULL;
  bmo_polar = NULL;
  bmo_atoms = NULL;
  bmo_num = NULL;
  vmo_occ = NULL;
  vmo_atoms = NULL;
  vmo_num = NULL;
  wAO = NULL;
  tAO = NULL;

  nb = 0;
  blist = NULL;

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz[i] = xyz_i[i];  
  for (int i=0;i<natoms;i++)
    q[i] = 0.;

  return;
}

void NBO::reset(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i)
{
  hasNBO = 0;
  natoms = natoms_i;
  bmo = 0;

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz[i] = xyz_i[i];    
  for (int i=0;i<natoms;i++)
    q[i] = 0.;

  return;
}

void NBO::freemem()
{
  hasNBO = 0;
  natoms = 0;
  bmo = 0;

  delete [] xyz;
  delete [] anumbers;
  delete [] anames;
  delete [] q;
  if (MO!=NULL)
    delete [] MO;
  if (mo_occ!=NULL)
    delete [] mo_occ;
  if (bmo_occ!=NULL)
    delete [] bmo_occ;
  if (bmo_atoms!=NULL)
    delete [] bmo_atoms;
  if (bmo_polar!=NULL)
    delete [] bmo_polar;
  if (bmo_num!=NULL)
    delete [] bmo_num;
  if (vmo_occ!=NULL)
    delete [] vmo_occ;
  if (vmo_atoms!=NULL)
    delete [] vmo_atoms;
  if (vmo_num!=NULL)
    delete [] vmo_num;

  return;
}



int NBO::read_mopac_mo(string filename)
{
  //printf("  reading MO from MOPAC file %s \n",filename.c_str());
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile){
    printf(" error opening MO file (MOPAC type): %s \n",filename.c_str());
    //exit(-1);
    return 0;
  }

 //build basis set
  nao = 0;
  mheadersize = 4 + natoms;
  for (int i=0;i<natoms;i++)
  {
    if (anumbers[i]==1)
    {
      nao += 1;
      mheadersize += 6;
    }
    else if (anumbers[i] > 1 && anumbers[i]<=10)
    {
      nao += 4;
      mheadersize += 6;
    }
    else
    {
      printf(" NBO: later row elements not supported yet! \n"); 
      exit(1);
    }
  }
  mheader = new string[mheadersize];
  mheader[0] = "[Molden Format]";
  mheader[1] = "[Atoms] (Angs)";
  char* str = new char[250];
  for (int i=0;i<natoms;i++)
  {
    sprintf(str," %2s %2i %2i %8.5f %8.5f %8.5f ",anames[i].c_str(),i+1,anumbers[i],xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]);
    mheader[i+2] = str;
  }

  int mc = 2 + natoms;
  mheader[mc++] = "[GTO]";
  for (int i=0;i<natoms;i++)
  {
    sprintf(str," %2i  0",i+1);
    mheader[mc++] = str;
    if (anumbers[i]==1)
    {
      mheader[mc++] = "S    3    1.000000";
      mheader[mc++] = "   3.42525091E+00    1.54328970E-01 ";
      mheader[mc++] = "   6.23913730E-01    5.35328140E-01 ";
      mheader[mc++] = "   1.68855400E-01    4.44634540E-01 ";
      mheader[mc++] = " ";
    }
    else if (anumbers[i]==5)
    {
      mheader[mc++] = "SP   3    1.000000";
      mheader[mc++] = "   2.23695610E+00   -9.99672300E-02   1.55916270E-01 ";
      mheader[mc++] = "   5.19820500E-01    3.99512830E-01   6.07683720E-01 ";
      mheader[mc++] = "   1.69061800E-01    7.00115470E-01   3.91957390E-01 ";
      mheader[mc++] = " ";
    }
    else if (anumbers[i]==6)
    {
     // mheader[mc++] = "S    3    1.000000";
     // mheader[mc++] = "   7.16168370E+01    1.54328970E-01 ";
     // mheader[mc++] = "   1.30450960E+01    5.35328140E-01 ";
     // mheader[mc++] = "   3.53051220E+00    4.44634540E-01 ";
      mheader[mc++] = "SP   3    1.000000";
      mheader[mc++] = "   2.94124940E+00   -9.99672300E-02   1.55916270E-01 ";
      mheader[mc++] = "   6.83483100E-01    3.99512830E-01   6.07683720E-01 ";
      mheader[mc++] = "   2.22289900E-01    7.00115470E-01   3.91957390E-01 ";
      mheader[mc++] = " ";
    }
    else if (anumbers[i]==7)
    {
     // mheader[mc++] = "S    3    1.000000";
     // mheader[mc++] = "   9.91061690E+01    1.54328970E-01 ";
     // mheader[mc++] = "   1.80523120E+01    5.35328140E-01 ";
     // mheader[mc++] = "   4.88566020E+00    4.44634540E-01 ";
      mheader[mc++] = "SP   3    1.000000";
      mheader[mc++] = "   3.78045590E+00   -9.99672300E-02   1.55916270E-01 ";
      mheader[mc++] = "   8.78496600E-01    3.99512830E-01   6.07683720E-01 ";
      mheader[mc++] = "   2.85714400E-01    7.00115470E-01   3.91957390E-01 ";
      mheader[mc++] = " ";
    }
    else if (anumbers[i]==8)
    {
      mheader[mc++] = "SP   3    1.000000";
      mheader[mc++] = "   5.03315130E+00   -9.99672300E-02   1.55916270E-01 ";
      mheader[mc++] = "   1.16959610E+00    3.99512830E-01   6.07683720E-01 ";
      mheader[mc++] = "   3.80389000E-01    7.00115470E-01   3.91957390E-01 ";
      mheader[mc++] = " ";
    }
    else if (anumbers[i]==9)
    {
      mheader[mc++] = "SP   3    1.000000";
      mheader[mc++] = "   6.46480320E+00   -9.99672300E-02   1.55916270E-01 ";
      mheader[mc++] = "   1.50228120E+00    3.99512830E-01   6.07683720E-01 ";
      mheader[mc++] = "   4.88588500E-01    7.00115470E-01   3.91957390E-01 ";
      mheader[mc++] = " ";
    }
    else 
    {
      printf("  element # %i not implemented! \n",anumbers[i]);
      exit(1);
    }
  }
  mheader[mc++] = "[MO]";


  string line;
  bool success=true;
  vector<string> tok_line;

 //allocate for all NBO and MO data
  alloc_nbo();
  MO = new double[nao*nao];
  for (int i=0;i<nao*nao;i++) MO[i] = 0.;

 //get MO's
  int found1 = 0;
  int cao = 0;
  int offset = 0;
  while(!infile.eof())
  {
    success = getline(infile,line);
    if (line.find("NO. OF FILLED LEVELS")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      nalpha = atoi(tok_line[5].c_str());
      nbeta = nalpha;
      nelec = nalpha + nbeta;
    }
    if (line.find("LOCALIZED ORBITALS")!=string::npos)
    {
      //cout << " RR: " << line << endl;
      found1 = 1;
    }
    if (found1)
    if (line.find("ROOT NO")!=string::npos)
    {
     // printf(" found ROOT NO, cao: %2i \n",cao);
      //cout << " RR1: " << line << endl;
      success = getline(infile,line);
      //cout << " RR2: " << line << endl;
      success = getline(infile,line);
      //cout << " RR3: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");

      //this line is orbital energies
      int nnew = tok_line.size();
      //tok_line = StringTools::tokenize(line, " \t");
      //printf(" about to write mo_occ. nnew: %2i offset: %2i \n",nnew,offset); fflush(stdout);
      for (int i=0;i<nnew;i++)
        mo_occ[offset+i] = atof(tok_line[i].c_str());

      success = getline(infile,line);
      //cout << " RR4: " << line << endl;
      success = getline(infile,line);
      //cout << " RR5: " << line << endl;

      //starting here, MO coefs
      for (int i=0;i<natoms;i++)
      {
        success = getline(infile,line);
        //cout << " RRat: " << line << endl;
        tok_line = StringTools::tokenize(line, " \t");
        nnew = tok_line.size() - 3;
        if (tok_line[1]=="H")
        {
         // printf(" found H \n");
         // printf(" offset: %2i cao: %2i nnew: %2i \n",offset,cao,nnew);
          tok_line = StringTools::tokenize(line, " \t");
          for (int j=0;j<nnew;j++)
            MO[(offset+j)*nao+cao] = atof(tok_line[3+j].c_str());
          cao += 1;
          success = getline(infile,line);
        }
        else if (tok_line[1]=="B" || tok_line[1]=="C" || tok_line[1]=="N" || tok_line[1]=="O" || tok_line[1]=="F")
        {
         // printf(" found C \n");
         // printf(" C offset: %2i cao: %2i nnew: %2i \n",offset,cao,nnew);
          for (int n=0;n<4;n++)
          {
            tok_line = StringTools::tokenize(line, " \t");
            for (int j=0;j<nnew;j++)
              MO[(offset+j)*nao+cao] = atof(tok_line[3+j].c_str());
            cao += 1;
            success = getline(infile,line);
          }
        }
      }
      cao = 0;
      offset += 6;
    } //if "ROOT NO" && found1
  }
  infile.close();

  nmo = nao;

  delete [] str;


#if 0
  printf("  printing MO coeff matrix (%2i mo's %2i ao's) \n",nmo,nao);
  for (int i=0;i<nmo;i++)
  { 
    for (int j=0;j<nao;j++)
      printf(" %8.5f",MO[i*nao+j]);
    printf("\n");
  }
#endif


 //now get other NBO data
  infile.open(filename.c_str());
  int cmo = 0;
  int nhp = 0;
  while (!infile.eof())
  {
    success = getline(infile,line);
    if (line.find("NET ATOMIC CHARGES")!=string::npos)
    {
      success = getline(infile,line);
      success = getline(infile,line);
      for (int i=0;i<natoms;i++)
      {
        success = getline(infile,line);
       // cout << " RR1: " << line << endl;
        tok_line = StringTools::tokenize(line, " \t");
        q[i] = atof(tok_line[2].c_str());
        double s1 = atof(tok_line[4].c_str());
        if (fabs(s1)<0.01) s1 = 0.01;
        double s2 = 0.;
        if (tok_line.size()>5)
          s2 = atof(tok_line[5].c_str());
        spratio[i] = s2/s1;
        //printf("  spratio[%2i]: %8.5f \n",i+1,spratio[i]);
      }
    }
    if (line.find("NUMBER OF CENTERS")!=string::npos)
    {
      success = getline(infile,line);
      success = getline(infile,line);
      while (!infile.eof())
      {
        success = getline(infile,line);
        tok_line = StringTools::tokenize(line, " \t");
        int lsize = tok_line.size();
        if (lsize>0)
        {
          if (line.find("LOCALIZED ORBITALS")!=string::npos)
            break;

          bmo_occ[cmo] = atof(tok_line[0].c_str());
          if (lsize==3 && bmo_occ[cmo]>100.) //was bmo_occ[cmo] == 1.0, handling bad mopac printing
          {
            tok_line = StringTools::tokenize(line, " .");
            bmo_atoms[2*cmo+0] = atoi(tok_line[2].c_str())/1000;
            bmo_polar[cmo] = -1.0;
            bmo_num[cmo] = cmo+1;
            //printf(" bmo_num: %i (cmo: %i) \n",bmo_num[cmo],cmo);
            //if (abs(bmo_num[bmo])>100) printf(" WARNING: bad bmo_num: %i \n",bmo_num[cmo]);
            cmo++;
          }
          else if (lsize==4)
          {
          // cout << " RR1: " << line << endl;
            bmo_atoms[2*cmo+0] = atoi(tok_line[2].c_str());
            bmo_polar[cmo] = -1.0;
            bmo_num[cmo] = cmo+1;
            //printf(" bmo_num: %i (cmo: %i) \n",bmo_num[cmo],cmo);
            //if (abs(bmo_num[bmo])>100) printf(" WARNING: bad bmo_num: %i \n",bmo_num[cmo]);
            cmo++;
          }
          else if (lsize==7)
          {
            bmo_atoms[2*cmo+0] = atoi(tok_line[2].c_str());
            bmo_atoms[2*cmo+1] = atoi(tok_line[5].c_str());
            double pl = atof(tok_line[3].c_str()) - atof(tok_line[6].c_str());
            bmo_polar[cmo] = -pl/100.;
            bmo_num[cmo] = cmo+1;
            //printf(" bmo_num: %i (cmo: %i) \n",bmo_num[cmo],cmo);
            //if (abs(bmo_num[cmo])>100) printf(" WARNING: bad bmo_num: %i \n",bmo_num[cmo]);

#if 0
           //for now, leave in highly polarized cases
            if (fabs(pl)>95.)
            {
              //printf("  highly polarized bond: %2i - %2i \n",bmo_atoms[2*cmo+0],bmo_atoms[2*cmo+1]);
              bmo_atoms[2*cmo+1] = 0;
              nhp++;
            }
#endif
            cmo++;
          }
          else if (lsize>9)
          {
            get_three_center(cmo,tok_line);
          }
        } //if line contains data
      } //while reading
    } //found localized orbital data

  } //while reading
  //printf("  found %2i highly polarized bonds \n",nhp);
  bmo = cmo;
 
  return 1;
}

void NBO::get_three_center(int& cmo, vector<string> tok_line)
{
  double wb3 = atof(tok_line[9].c_str());
 //third contributor is small
  if (wb3<10.)
  {
    bmo_atoms[2*cmo+0] = atoi(tok_line[2].c_str());
    bmo_atoms[2*cmo+1] = atoi(tok_line[5].c_str());
    double pl = atof(tok_line[3].c_str()) - atof(tok_line[6].c_str());
    bmo_polar[cmo] = -pl/100.;
    bmo_num[cmo] = cmo+1;
    cmo++;
    return;
  }
  else
  {
    int b1,b2,b3;
    b1 = atoi(tok_line[2].c_str());
    b2 = atoi(tok_line[5].c_str());
    b3 = atoi(tok_line[8].c_str());
    double p1a = atof(tok_line[3].c_str());
    double p2a = atof(tok_line[6].c_str());
    double p3a = atof(tok_line[9].c_str());
    double pl1 = p1a - p2a;
    double pl2 = p2a - p3a;

   //put H in the middle
    if (anumbers[b2-1]!=1)
    {
      int t1 = b2;
      if (anumbers[b1-1]==1)
      {
        b2 = b1;
        b1 = t1;
        pl1 = p2a - p1a;
        pl2 = p1a - p3a;
      }
      if (anumbers[b3-1]==1)
      {
        b2 = b3;
        b3 = t1;
        pl1 = p1a - p3a;
        pl2 = p3a - p2a;
      }
    }
    bmo_atoms[2*cmo+0] = b1;
    bmo_atoms[2*cmo+1] = b2;
    bmo_num[cmo] = cmo+1;
    bmo_polar[cmo] = -pl1/100.;
    cmo++;
    bmo_atoms[2*cmo+0] = b2;
    bmo_atoms[2*cmo+1] = b3;
    bmo_num[cmo] = cmo+1;
    bmo_polar[cmo] = -pl2/100.;
    cmo++;
  }

  return;
}

int NBO::read_mo(string filename)
{
  printf("  reading MO from Molden file %s. ",filename.c_str());
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile){
    printf(" error opening MO file (Molden formatted): %s \n",filename.c_str());
    return 0;
  }

#if ISOMER_DB || 1
  alloc_nbo();
  printf(" read_mo not working right now, not using \n");
  return 1;
#endif

  string line;
  bool success=true;    
  vector<string> tok_line;

  int MAX_BASIS = 25 * natoms;
  wAO = new int[MAX_BASIS];
  tAO = new int[MAX_BASIS];
  int w = 0;
  nao = 0;

 //save header information
  int mcounting = 0;
  mheadersize = 0;
  while(!infile.eof())
  {
    success = getline(infile,line);
    if (line.find("[MO]")!=string::npos)
      break;
    if (mcounting)
      mheadersize++;
    if (line.find("[Molden Format]")!=string::npos)
      mcounting = 1;
  }
  infile.close();

  int mc = 0;
  mcounting = 0;
  mheader = new string[mheadersize+2];
  infile.open(filename.c_str());
  while(!infile.eof())
  {
    success = getline(infile,line);
    if (line.find("[Molden Format]")!=string::npos)
      mcounting = 1;
    if (mcounting)
      mheader[mc++] = line;
    if (mc>=mheadersize) break;
  }
  mheader[mc++] = " ";
  mheader[mc++] = "[MO]";
  mheadersize += 2;
  infile.close();

 //get basis and MO data
  infile.open(filename.c_str());
  while(!infile.eof())
  {
    success = getline(infile,line);
    if (line.find("[GTO]")!=string::npos)
    {
     // cout << " RR00: " << line << endl;
      while(line.find("[MO]")==string::npos)
      {
        success = getline(infile,line);
       // cout << " RR0: " << line << endl;
        tok_line = StringTools::tokenize(line, " \t");
        if (tok_line.size()>1)
        {
          int isAO = atoi(tok_line[1].c_str());
          if (isAO==0)
          {
            int cont = 1;
            while(cont)
            {
              success = getline(infile,line);
              tok_line = StringTools::tokenize(line, " \t");
              //cout << "  RR0: " << line << endl;
              if (line.find("S")!=string::npos)
              {
                tAO[nao] = 0;
                wAO[nao++] = w;
              }
              if (line.find("P")!=string::npos)
              {
                for (int i=0;i<3;i++)
                {
                  tAO[nao] = 1;
                  wAO[nao++] = w;
                }
              }
              if (line.find("D")!=string::npos) //spd only
              {
                for (int i=0;i<6;i++) //hardcoded: 6D functions
	        {
                  tAO[nao] = 2;
                  wAO[nao++] = w;
                }
              }
              if (tok_line.size()<1)
                break;
            }
            w++;
          }
        }
      } //Reading AO basis
      if (nao>MAX_BASIS) 
      {
        printf(" nao exceeds MAX_BASIS in nbo.cpp! nao: %3i MAX_BASIS: %3i \n",nao,MAX_BASIS);
        exit(1);
      }
      break;
    }
  }
  for (int i=0;i<4;i++)
    success = getline(infile,line);

  printf("  found %i ao's (with 3p 5d etc) \n",nao);

  double* jMO1 = new double[MAX_BASIS]; //Note max basis limit
  int dim = 0;
  while (!infile.eof())
  {
    success = getline(infile,line);
    if (line.find("Sym")!=string::npos)
    {
      //cout << " RR1: " << line << endl;
      break;
    }
    tok_line = StringTools::tokenize(line, " \t");
    jMO1[dim] = atof(tok_line[1].c_str());
    dim++;
    if (dim>nao)
    {
      printf("  not expecting this many orbitals: %2i \n",dim);
      exit(1);
    }
  }
  printf(" dimension: %3i \n",dim);
  nao = dim;

  MO = new double[dim*nao];
  for (int i=0;i<dim*nao;i++) MO[i] = 0.;
  for (int i=0;i<nao;i++) MO[i] = jMO1[i];

  nmo = 1;
  while (!infile.eof())
  {
    success = getline(infile,line);
    cout << " RRX: " << line << endl;
    if (line.find("Occup=")!=string::npos)
    {
      for (int j=0;j<nao;j++)
      {
        success = getline(infile,line);
        cout << " RRY: " << line << endl;
        tok_line = StringTools::tokenize(line, " \t");
        MO[nmo*nao+j] = atof(tok_line[1].c_str());
      } //loop over c MO coefs
      nmo++;
    } //if found MO
  }

  infile.close();

#if 0
  printf("  printing AO-atom matching \n");
  for (int i=0;i<nao;i++)
    printf(" %2i",wAO[i]);
  printf("\n");
  printf("  printing AO-type \n");
  for (int i=0;i<nao;i++)
  {
    if (tAO[i]==0)
      printf("  s");
    else if (tAO[i]==1)
      printf("  p");
    else if (tAO[i]==2)
      printf("  d");
  }
  printf("\n");
#endif
#if 0
  printf("  printing MO coeff matrix (%i mo's) \n",nmo);
  for (int i=0;i<nmo;i++)
  { 
    for (int j=0;j<nao;j++)
      printf(" %8.5f",MO[i*nao+j]);
    printf("\n");
  }
#endif

  delete [] jMO1;

  return dim;
}
