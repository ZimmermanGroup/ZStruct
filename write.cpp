// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "zstruct.h"


int ZStruct::get_frozen(int nreact1)
{
  int natoms_max = 200; //same in update ftn
  update_active_mem(nreact1);

  ifstream frzfile;
  int nfound = 0;
  vector <string> tok_line;
  string line;
  bool success;
  struct stat sts;
  for (int i=1;i<=nreact1;i++)
  {
    string nstr = StringTools::int2str(i,1,"0");
    string filename = "frozen"+nstr+".xyz";
    if (stat(filename.c_str(), &sts) != -1)
    {
      //printf(" found file %s \n",filename.c_str());

      frzfile.open(filename.c_str());
      if (!frzfile)
      {
        printf(" ERROR: couldn't open: %s \n",filename.c_str());
        exit(1);
      }
      
      int nfound = 0;
      while (!frzfile.eof())
      {
        success=getline(frzfile, line);
        tok_line = StringTools::tokenize(line, " \t");

        //printf(" RR: %s \n",line.c_str()); fflush(stdout);
 
        if (tok_line.size()==1)
        {
          int r1 = atoi(tok_line[0].c_str());
          int a1 = abs(r1) -1;
          if (r1>0)
            active[i-1][a1] = 0;
          else
            active[i-1][a1] = -1;
        }
        else if (tok_line.size()>1)
          printf("\n WARNING: each line in frozen.xyz must have 1 entry \n");     

      } //while reading frzfile

      frzfile.close();

      nfound++;
    }
  }

  printf("\n Inactive atoms: \n");
  for (int i=0;i<nreact1;i++)
  {
    printf(" reactant %2i:",i);
    for (int j=0;j<natoms_max;j++)
    if (active[i][j]<1)
      printf(" %2i",j+1);
    printf("\n");
  }
  printf("\n");

#if 0
  string frzfilestr = "frozen.xyz"; //later change to species basis
  frzfile.open(frzfilestr.c_str());
  if (!frzfile)
  {
    printf("\n WARNING: couldn't find frozen.xyz file \n");
    return 0;
  }

  vector <string> tok_line;
  string line;
  bool success;
  int nfound = 0;
  while (!frzfile.eof())
  {
    success=getline(frzfile, line);
    tok_line = StringTools::tokenize(line, " \t");

    //printf(" RR: %s \n",line.c_str()); fflush(stdout);
 
    if (tok_line.size()==1)
    {
      int num = atoi(tok_line[0].c_str()) -1;
      active[num] = 0;
      nfound++;
    }
    else if (tok_line.size()>1)
      printf("\n WARNING: each line in frozen.xyz must have 1 entry \n");
     
  }
  printf("\n Inactive atoms:");
  for (int i=0;i<natoms_max;i++)
  if (active[i]==0)
  {
    printf(" %2i",i+1);
  }
  printf("\n");
#endif

  return nfound;
}


void ZStruct::update_active_mem(int nreact1)
{
  if (nreact1<activea) return;

  int natoms_max = 200;
  if (active==NULL)
  {  
    active = new int*[nreact1];
    for (int i=0;i<nreact1;i++)
      active[i] = new int[natoms_max];
    for (int i=0;i<nreact1;i++)
    for (int j=0;j<natoms_max;j++)
      active[i][j] = 1;
  }
  else
  {
    int** active1 = new int*[nreact1];
    for (int i=0;i<activea;i++)
      active1[i] = active[i];
    for (int i=activea;i<nreact1;i++)
      active1[i] = new int[natoms_max];
    for (int i=activea;i<nreact1;i++)
    for (int j=0;j<natoms_max;j++)
      active1[i][j] = 1;

    delete [] active;
    active = active1;  
  }
  activea = nreact1;

  return;
}

void ZStruct::assign_sequential_pairs(int wp)  
{
 //deprecated
  //printf(" in assign_sequential_pairs wp: %i \n",wp);
  
  npair = wp+1;
          
  int nfound = 0;
  for (int i=0;i<nreact;i++)
  for (int j=0;j<=i;j++)
  {
    if (nfound>=npair)
      break;
    pairs[2*nfound+0] = i;
    pairs[2*nfound+1] = j;
    nfound++;
  }       

  if (nreact==1)
  {
    pairs[2*wp+0] = 0;
    pairs[2*wp+1] = 0;
  }

  for (int i=0;i<npair;i++)
    printf(" pairs(asp): %i %i \n",pairs[2*i+0],pairs[2*i+1]);

  return;
}

void ZStruct::save_shuttles(int type)
{
#if DO_NOT_WRITE && 1
  printf("\n not saving GSMDATA_SH file \n");
  return;
#endif

  if (nshuttle<1) return;

  string sh_string = "GSMDATA_SH";
  if (type==0)
    sh_string = "GSMDATA_SH_sub";

  ofstream shfile;
  shfile.open(sh_string.c_str());

  for (int i=0;i<niso;i++)
  {
    shfile << " " << wshpair[i] << endl;
  }

  shfile.close();

  return;
}


void ZStruct::save_gsm(int type)
{
#if DO_NOT_WRITE
  printf("  not saving GSMDATA file \n");
  printf("  would have written: \n");

  for (int i=0;i<nreact;i++)
  {
    if (reactnum[2*i+0]!=reactnum[2*i+1])
    {
      cout << " REACT " << i << "   " 
              << reactnum[2*i+0] << " " << reactnum[2*i+1]-1 << endl;
    }
    else
    {
      cout << " REACT " << i << "  " 
              << reactnum[2*i+0] << " -1" << endl;
    }
  }

  for (int i=0;i<npair;i++)
  {
    if (pairsnum[2*i+0]!=pairsnum[2*i+1])
    {
      cout << " PAIR " << pairs[2*i+0] << " " << pairs[2*i+1] << "  " 
              << pairsnum[2*i+0] << " " << pairsnum[2*i+1]-1 << endl;
    }
    else
    {
      cout << " PAIR " << pairs[2*i+0] << " " << pairs[2*i+1] << "  " 
              << pairsnum[2*i+0] << " -1" << endl;
    }
  }

  return;
#endif

//  printf(" save_gsm, nreact: %i npair: %i \n",nreact,npair);

  string gsm_string = "GSMDATA";
  if (type==0)
    gsm_string = "GSMDATA_sub";

  ofstream gsmfile;
  gsmfile.open(gsm_string.c_str());

  for (int i=0;i<nreact;i++)
  {
    if (reactnum[2*i+0]!=reactnum[2*i+1])
    {
      gsmfile << " REACT " << i << "   " 
              << reactnum[2*i+0] << " " << reactnum[2*i+1]-1 << endl;
    }
    else
    {
      gsmfile << " REACT " << i << "  " 
              << reactnum[2*i+0] << " -1" << endl;
    }
  }

  for (int i=0;i<npair;i++)
  {
    if (pairsnum[2*i+0]!=pairsnum[2*i+1])
    {
      gsmfile << " PAIR " << pairs[2*i+0] << " " << pairs[2*i+1] << "  " 
              << pairsnum[2*i+0] << " " << pairsnum[2*i+1]-1 << endl;
    }
    else
    {
      gsmfile << " PAIR " << pairs[2*i+0] << " " << pairs[2*i+1] << "  " 
              << pairsnum[2*i+0] << " -1" << endl;
    }
  }
#if 0
  for (int i=0;i<nsave;i++)
  {
    gsmfile << "START " << start[i] << endl;
    gsmfile << "END " << end[i] << endl;
  }
#endif

  gsmfile.close();

  return;
}

void ZStruct::save_frozen(int natoms1, int* active1, int wfile)
{
#if SKIPGSM || DO_NOT_WRITE
  printf(" not saving frozen file \n");
  return;
#endif

  //printf(" save_frozen, natoms1: %i wfile: %i \n",natoms1,wfile);

  string nstr = StringTools::int2str(wfile,1,"0");
  string frzstr = "frozen"+nstr+".xyz";

  ofstream frzfile;
  frzfile.open(frzstr.c_str());

  for (int i=0;i<natoms1;i++)
  if (!active1[i])
  {
    frzfile << " " << i+1 << endl;
  }

  frzfile.close();

  return;
}



void ZStruct::write_initial_xyz(int wfile, int natoms1, string* anames1, double* xyz1, int q1)
{
#if DO_NOT_WRITE
  return;
#endif

  //printf(" in write_initial_xyz \n"); fflush(stdout);

  string id = StringTools::int2str(wfile,4,"0");  
  string xyzfile_string = "scratch/initial"+id+".xyz";

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(8);
  xyzfile << fixed;

  int natoms2 = natoms1;
  for (int i=0;i<natoms2;i++)
  if (anames1[i]=="X")
    natoms2--;

  xyzfile << natoms2 << endl;
  xyzfile << " " << q1 << endl;
  for (int i=0;i<natoms1;i++)
  if (anames1[i]!="X")
    xyzfile << " " << anames1[i] << " " << xyz1[3*i+0] << " " << xyz1[3*i+1] << " " << xyz1[3*i+2] << endl;

  xyzfile.close();

  return;
}

void ZStruct::write_pair_xyz(int wfile, int natoms1, string* anames1, double* xyz1)
{
#if DO_NOT_WRITE
  return;
#endif

  //printf(" in write_pair_xyz \n"); fflush(stdout);

  string id = StringTools::int2str(wfile,4,"0");  
  string xyzfile_string = "scratch/rpair"+id+".xyz";

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(8);
  xyzfile << fixed;

  int natoms2 = natoms1;
  for (int i=0;i<natoms2;i++)
  if (anames1[i]=="X")
    natoms2--;

  xyzfile << natoms2 << endl << endl;
  for (int i=0;i<natoms1;i++)
  if (anames1[i]!="X")
    xyzfile << " " << anames1[i] << " " << xyz1[3*i+0] << " " << xyz1[3*i+1] << " " << xyz1[3*i+2] << endl;

  xyzfile.close();

  return;
}

void ZStruct::write_ISOMER(int wfile, int nadd, int* add, int nbrk, int* brk)
{
  //printf(" writing ISOMER %2i \n",wfile);
  write_ISOMER(wfile,nadd,add,nbrk,brk,0,NULL,NULL,0,NULL,NULL);
}

void ZStruct::write_ISOMER(int wfile, int nadd, int* add, int nbrk, int* brk, int nangle, int** angles, double* anglev, int ntor, int** torsion, double* torv)
{
#if DO_NOT_WRITE
  return;
#endif

  //printf(" in write_ISOMER \n"); fflush(stdout);

  string id = StringTools::int2str(wfile,4,"0");
  string ifile_string = "scratch/ISOMERS"+id;

  ofstream ifile;
  ifile.open(ifile_string.c_str());
  ifile << setprecision(8);

  ifile << "NEW" << endl;

  for (int i=0;i<nfragb;i++)
    ifile << " BOND " << fragb[2*i+0]+1 << " " << fragb[2*i+1]+1 << endl;

  for (int i=0;i<nadd;i++)
    ifile <<  " ADD "  << add[2*i+0]+1 << " " << add[2*i+1]+1 << endl;
  for (int i=0;i<nbrk;i++)
    ifile << " BREAK " << brk[2*i+0]+1 << " " << brk[2*i+1]+1 << endl;

  for (int i=0;i<nangle;i++)
    ifile << " ANGLE " << angles[i][0]+1 << " " << angles[i][1]+1 << " " << angles[i][2]+1 << " " << anglev[i] << endl;
  for (int i=0;i<ntor;i++)
    ifile << " TORSION " << torsion[i][0]+1 << " " << torsion[i][1]+1 << " " << torsion[i][2]+1 << " " << torsion[i][3]+1 << torv[i] << endl;

  ifile << endl;

  ifile.close();


#if 0
  cout << endl << ifile_string << endl;
  cout << "NEW" << endl;

  for (int i=0;i<nadd;i++)
    cout <<  " ADD "  << add[2*i+0]+1 << " " << add[2*i+1]+1 << endl;
  for (int i=0;i<nbrk;i++)
    cout << " BREAK " << brk[2*i+0]+1 << " " << brk[2*i+1]+1 << endl;

  for (int i=0;i<nangle;i++)
    cout << " ANGLE " << angles[3*i+0]+1 << " " << angles[3*i+1]+1 << " " << angles[3*i+2]+1 << " " << anglev[i]+1 << endl;
  for (int i=0;i<ntor;i++)
    cout << " TORSION " << torsion[4*i+0]+1 << " " << torsion[4*i+1]+1 << " " << torsion[4*i+2]+1 << " " << torsion[4*i+3]+1 << torv[i] << endl;

  cout << endl;
#endif

  return;
}


int ZStruct::read_ISOMER(int wfile, int& nadd, int* add, int& nbrks, int* brks)
{
  nadd = nbrks = 0;

  string nstr=StringTools::int2str(wfile,4,"0");
  string isomerfile = "scratch/ISOMERS"+nstr;

  ifstream output(isomerfile.c_str(),ios::in);
  if (!output)
  {
    isomerfile = "scratch/inputs/ISOMERS"+nstr;
    output.open(isomerfile.c_str());
    if (!output)
    {
      printf(" couldn't find ISOMERS file \n");
      return 0;
    }
  }

  string line;
  vector<string> tok_line;
  while(!output.eof())
  {
    getline(output,line);
    //cout << " RR " << line << endl;

    if (line.find("ADD")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      add[2*nadd] = atoi(tok_line[1].c_str()) -1;
      add[2*nadd+1] = atoi(tok_line[2].c_str()) -1;
      //printf(" add: %i %i \n",add[2*nadd]+1,add[2*nadd+1]+1);
      nadd++;
    }
    if (line.find("BREAK")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      brks[2*nbrks] = atoi(tok_line[1].c_str()) -1;
      brks[2*nbrks+1] = atoi(tok_line[2].c_str()) -1;
      //printf(" brk: %i %i \n",brks[2*nbrks]+1,brks[2*nbrks+1]+1);
      nbrks++;
    }
#if 0
    if (line.find("ANGLE")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      angles[3*nangle+0] = atoi(tok_line[1].c_str()) -1;
      angles[3*nangle+1] = atoi(tok_line[2].c_str()) -1;
      angles[3*nangle+2] = atoi(tok_line[3].c_str()) -1;
      anglet[nangle] = atof(tok_line[4].c_str());
      printf(" angle: %i %i %i align to %4.3f \n",angles[3*nangle+0]+1,angles[3*nangle+1]+1,angles[3*nangle+2]+1,anglet[nangle]);
      nangle++;
    }
    if (line.find("TORSION")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      tors[4*ntors+0] = atoi(tok_line[1].c_str()) -1;
      tors[4*ntors+1] = atoi(tok_line[2].c_str()) -1;
      tors[4*ntors+2] = atoi(tok_line[3].c_str()) -1;
      tors[4*ntors+3] = atoi(tok_line[4].c_str()) -1;
      tort[ntors] = atof(tok_line[5].c_str());
      printf(" tor: %i %i %i %i align to %4.3f \n",tors[4*ntors+0]+1,tors[4*ntors+1]+1,tors[4*ntors+2]+1,tors[4*ntors+3]+1,tort[ntors]);
      ntors++;
    }
#endif
  } //while !eof

#if 0
  printf(" ISOMER%s \n",nstr.c_str());
  for (int i=0;i<nadd;i++)
    printf("   adding: %2i %2i \n",add[2*i+0]+1,add[2*i+1]+1);
  for (int i=0;i<nbrks;i++)
    printf("   breaking: %2i %2i \n",brks[2*i+0]+1,brks[2*i+1]+1);
#endif

  int success = 1;
  if (nadd+nbrks<1) success = 0;

  return success;
}


int ZStruct::read_reactants(ICoord* icr)
{
  int nfound = 0;

  int maxr = 500;
  struct stat sts;

  int* has_tm = new int[maxr];
  for (int i=0;i<maxr;i++) has_tm[i] = 0;

  for (int i=1;i<=maxr;i++)
  {
    string nstr = StringTools::int2str(i,1,"0");
    string filename = "react"+nstr+".xyz";
    if (stat(filename.c_str(), &sts) != -1)
    {
      printf(" found file %s, saved to icr[%i] \n",filename.c_str(),nfound);
      icr[nfound].init(filename);
      icr[nfound].print_bonds();
      vector<string> tok_line = StringTools::tokenize(icr[nfound].comment, " \t");
      if (tok_line.size()>0)
      {
        //printf("  comment line first element: %s \n",tok_line[0].c_str());
        icr[nfound].q1 = atoi(tok_line[0].c_str());
        if (tok_line.size()>1)
          icr[nfound].s1 = atoi(tok_line[1].c_str());
      }
      for (int j=0;j<icr[nfound].natoms;j++)
      if (icr[nfound].isTM(j))
        has_tm[nfound]++;
      nfound++;
    }
    else
      break;
  }
  if (nfound==0)
  {
    printf(" didn't find react.xyz files, looking for test.xyz \n");
    string filename = "test.xyz";
    if (stat(filename.c_str(), &sts) != -1)
    {
      icr[0].init(filename);
      icr[0].print_bonds();
      icr[0].make_frags();
      create_frag_bonds(icr[0]);
      for (int j=0;j<icr[0].natoms;j++)
      if (icr[0].isTM(j))
        has_tm[0]++;
      nfound++;
    }
  }

  int nreact0 = nreact;
  nreact = nfound;
  update_react_mem();
  nreact = nreact0;

  for (int i=0;i<nfound;i++) 
  {
    natomsr[i] = icr[i].natoms;
    if (has_tm[i])
      pair_react[i] = 0;
  }
#if TM_PAIRS
  printf(" allowing metal pairs! \n");
  for (int i=0;i<nfound;i++)
    pair_react[i] = 1;
#endif

  delete [] has_tm;

  return nfound;
}


int ZStruct::read_string(int wfile, double* energies, int natoms0, double** xyz)
{
  int ngeom = 0;

  string nstr=StringTools::int2str(wfile,4,"0");
  string strfilename = "stringfile.xyz"+nstr;

  ifstream strfile;
  strfile.open(strfilename.c_str());
  if (!strfile)
  {
    strfilename = "savestrings/stringfile.xyz"+nstr;
    strfile.open(strfilename.c_str());
    if (!strfile)
    {
      //printf(" couldn't find %s file \n",strfilename.c_str());
      return 0;
    }
  }

  string line;
  bool success=true;
  success=getline(strfile, line);
  int length=StringTools::cleanstring(line);
  int natoms1 = atoi(line.c_str());

  if (natoms1!=natoms0)
  {
    printf(" natom mismatch in read_string found: %i expected: %i file: %s \n",natoms1,natoms0,strfilename.c_str());
    return -natoms1;
//    exit(1);
  }

  int first = 1;
  while(!strfile.eof())
  {
    if (!first) success = getline(strfile, line);
    else first = 0;

    success = getline(strfile, line);
    if (strfile.eof())
      break;

    int length = StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    energies[ngeom] = atof(tok_line[0].c_str());

    if (xyz!=NULL)
    for (int j=0;j<natoms1;j++)
    {
      success = getline(strfile, line);
      length = StringTools::cleanstring(line);
      tok_line = StringTools::tokenize(line, " \t");
//      cout << " i: " << i << " string: " << line << endl;
      xyz[ngeom][3*j+0] = atof(tok_line[1].c_str());
      xyz[ngeom][3*j+1] = atof(tok_line[2].c_str());
      xyz[ngeom][3*j+2] = atof(tok_line[3].c_str());
      if (strfile.eof()) { ngeom--; break; }
    }
    if (xyz==NULL)
    for (int j=0;j<natoms1;j++)
      success = getline(strfile, line);

    ngeom++;
  } //while !eof

  strfile.close();

  int nanfound = 0;
  if (xyz!=NULL)
  for (int i=0;i<ngeom;i++)
  if (!nanfound)
  for (int j=0;j<3*natoms1;j++)
  if (xyz[i][j]!=xyz[i][j])
  {
    nanfound = 1;
    break;
  }
  if (nanfound)
    return 0;

#if 0
  printf(" printing string back, file: %s \n",strfilename.c_str());
  if (xyz!=NULL)
  for (int i=0;i<ngeom;i++)
  {
    printf(" %i \n\n",natoms1);
    for (int j=0;j<natoms1;j++)
      printf(" %4.4f %4.4f %4.4f \n",xyz[i][3*j+0],xyz[i][3*j+1],xyz[i][3*j+2]);
  }  
#endif


  return ngeom;
}

double ZStruct::read_temperature()
{
  string filename = "TEMPERATURE";

  ifstream file;
  file.open(filename.c_str());
  if (!file)
  {
    printf(" couldn't find %s file \n",filename.c_str());
    return -1.;
  }

  string line;
  getline(file, line);
  int length=StringTools::cleanstring(line);
  double T = atof(line.c_str());

  file.close();

//  printf(" found T: %4.2f \n",T);

  return T;
}


#if 0
    if (line.find("START")!=string::npos)
    { 
      tok_line = StringTools::tokenize(line, " \t");
      start = atoi(tok_line[1].c_str());
    }
    else if (line.find("END")!=string::npos)
    { 
      tok_line = StringTools::tokenize(line, " \t");
      end = atoi(tok_line[1].c_str());
    }
#endif


int ZStruct::get_limits()
{
  nelem = 80;
  if (climit_l==NULL)
    climit_l = new int[nelem];
  for (int i=0;i<nelem;i++) climit_l[i] = 0;
  if (climit_h==NULL)
    climit_h = new int[nelem];
  for (int i=0;i<nelem;i++) climit_h[i] = -1;

  string limfilestr = "LIMITS";
  ifstream limfile;
  limfile.open(limfilestr.c_str());
  if (!limfile)
  {
    printf("\n ERROR: couldn't find LIMITS file \n");
    exit(-1);
  }

  vector <string> tok_line;
  string line;
  bool success;
  int nfound = 0;
  int anum;
  while (!limfile.eof())
  {
    success=getline(limfile, line);
    tok_line = StringTools::tokenize(line, " \t");

    //printf(" RR: %s \n",line.c_str()); fflush(stdout);
 
    if (line.find("ATOMS")!=string::npos)
    {
      atomlimit = atoi(tok_line[1].c_str());
      printf(" Found ATOMS limit of %2i \n",atomlimit);
    }
    else if (tok_line.size()==3)
    {
      anum = PTable::atom_number(tok_line[0]);
      int low = atoi(tok_line[1].c_str());
      int high = atoi(tok_line[2].c_str());
      if (low>high) 
      {
        int tmp = high;
        high = low;
        low = tmp;
      }
      climit_l[anum] = low;
      climit_h[anum] = high;
      nfound++;
    }
    else if (tok_line.size()!=0)
      printf("\n WARNING: each atom line in LIMITS must have 3 entries \n");
     
  }

  limfile.close();

#if 1
  printf("\n Coordination number limits \n");
  for (int i=1;i<nelem;i++)
  if (climit_h[i]>-1)
  {
    string element = PTable::atom_name(i);
    printf(" element %2s: %2i %2i \n",element.c_str(),climit_l[i],climit_h[i]);
  }
#endif

  rxns1.set_limits(climit_l,climit_h);
  rxns2.set_limits(climit_l,climit_h);

  return nfound;
}

int ZStruct::read_shuttles()
{
  int nfound = 0;
  int nmax = 10;
  if (wshuttle==NULL)
    wshuttle = new int[nmax];
  for (int i=0;i<nmax;i++) wshuttle[i] = -1;
  if (shuttles==NULL)
    shuttles = new int*[nmax];
  for (int i=0;i<nmax;i++)
    shuttles[i] = new int[4];

  string shfilestr = "SHUTTLES";
  ifstream shfile;
  shfile.open(shfilestr.c_str());
  if (!shfile)
  {
    printf("\n WARNING: couldn't find SHUTTLES file \n");
    return 0;
  }

  vector <string> tok_line;
  string line;
  bool success;
  int anum;
  while (!shfile.eof())
  {
    success=getline(shfile, line);
    tok_line = StringTools::tokenize(line, " \t");

    //printf(" RR: %s \n",line.c_str()); fflush(stdout);
 
    if (line.find("SHUTTLE")!=string::npos && tok_line.size()==6)
    {
      int a0 = atoi(tok_line[1].c_str())-1;
      int a1 = atoi(tok_line[2].c_str())-1;
      int a2 = atoi(tok_line[3].c_str())-1;
      int a3 = atoi(tok_line[4].c_str())-1;
      string e1 = tok_line[5];
      int a1e = PTable::atom_number(e1);
      printf(" found shuttle %i. H-O: %i-%i O: %i element: %s (%i) \n",a0+1,a1+1,a2+1,a3+1,e1.c_str(),a1e);
      wshuttle[nfound] = a0;
      shuttles[nfound][0] = a1;
      shuttles[nfound][1] = a2;
      shuttles[nfound][2] = a3;
      shuttles[nfound][3] = a1e;
      nfound++;
    }
    else if (tok_line.size()!=0)
    {
      printf("\n WARNING: each atom line in SHUTTLES must have X entries \n");
      printf(" Format: SHUTTLE  REACTNUM  H-O O  ELEMENT \n");
    }
    if (nfound>nmax)
    {
      printf(" ERROR: static memory alloc exceeded in read_shuttles: %i/%i \n",nfound,nmax);
      exit(1);
    }     
  }

  shfile.close();

  if (nfound>1)
  {
    printf(" ERROR: only one shuttle currently supported \n");
    exit(1);
  }

  return nfound;
}

int ZStruct::read_ga()
{
  string file = "scratch/ga.save";
  return rxns2.read_ga(file);
}
  
void ZStruct::write_ga()
{
  string file = "scratch/ga.save";
  rxns2.write_ga(file);
  return;
}


