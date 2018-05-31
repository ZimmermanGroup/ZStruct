// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "dft.h"
#include "constants.h"
//#include "utils.h"
using namespace std;

#define DO_NBO 1

// for para execution
void DFT::sp_dnr(string filename) {

  //printf(" in DFT/sp_dnr() \n");

  energy0 = energy = 0;

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  inpfile << " $molecule " << endl;
  inpfile << " " << DFT_CHARGE << " " << DFT_SPIN << endl;

  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  inpfile << " $end " << endl;
  inpfile << endl << endl;
  inpfile << " $rem " << endl;
  inpfile << " JOBTYPE SP " << endl;
  inpfile << " EXCHANGE " << DFT_EXCHANGE << endl;
  inpfile << " CORRELATION " << DFT_CORRELATION << endl;
#if UNRESTRICTED
  inpfile << " UNRESTRICTED TRUE " << endl;
#endif
  inpfile << " SCF_ALGORITHM rca_diis " << endl;
//  inpfile << " SCF_ALGORITHM gdm " << endl;
  inpfile << " SCF_MAX_CYCLES 250 " << endl;
  inpfile << " SCF_CONVERGENCE 6 " << endl;
#if B631GSS
  inpfile << " BASIS 6-31G** " << endl;
#elif B631GS
  inpfile << " BASIS 6-31G* " << endl;
#elif LANL2DZ
  inpfile << " BASIS LANL2DZ " << endl;
  inpfile << " ECP LANL2DZ " << endl;
#elif BASISMIX
  inpfile << " BASIS mixed " << endl;
  inpfile << " ECP gen " << endl;
#else
  inpfile << " BASIS 6-31G " << endl;     
#endif
//  inpfile << " WAVEFUNCTION_ANALYSIS FALSE " << endl;
#if DO_NBO
  inpfile << " NBO 1 " << endl;
  inpfile << " print_orbitals 250 " << endl;
  inpfile << " MOLDEN_FORMAT TRUE " << endl;
#endif

  inpfile << " $end " << endl;
  inpfile << endl;

#if (BASISMIX || BASISGEN)
  inpfile << endl;
  string cmd = "cat qmix >> "+inpfile_string;
  system(cmd.c_str());
  inpfile << endl;
#endif

//  string cmd = "qchem "+filename+" "+filename+".out";
//  system(cmd.c_str());

  return;

}

// for para execution
void DFT::opt_dnr(string filename) {

  //printf(" in DFT/opt_dnr() \n");

  energy0 = energy = 0;

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  inpfile << " $molecule " << endl;
  inpfile << " " << DFT_CHARGE << " " << DFT_SPIN << endl;

  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  inpfile << " $end " << endl;
  inpfile << endl << endl;
  inpfile << " $rem " << endl;
  inpfile << " JOBTYPE OPT " << endl;
  inpfile << " EXCHANGE " << DFT_EXCHANGE << endl;
  inpfile << " CORRELATION " << DFT_CORRELATION << endl;
#if UNRESTRICTED
  inpfile << " UNRESTRICTED TRUE " << endl;
#endif

  inpfile << " GEOM_OPT_MAX_CYCLES  " << DFT_OPT_CYCLES << endl;
  inpfile << " GEOM_OPT_TOL_DISPLACEMENT 2500 " << endl;
  inpfile << " GEOM_OPT_TOL_GRADIENT     800 " << endl;
  inpfile << " GEOM_OPT_TOL_ENERGY	5000 " << endl;

  inpfile << " SCF_ALGORITHM rca_diis " << endl;
//  inpfile << " SCF_ALGORITHM gdm " << endl;
  inpfile << " SCF_MAX_CYCLES 250 " << endl;
  inpfile << " SCF_CONVERGENCE 6 " << endl;
#if B631GSS
  inpfile << " BASIS 6-31G** " << endl;
#elif B631GS 
  inpfile << " BASIS 6-31G* " << endl;
#elif LANL2DZ
  inpfile << " BASIS LANL2DZ " << endl;
  inpfile << " ECP LANL2DZ " << endl;
#elif BASISMIX
  inpfile << " BASIS mixed " << endl;
  inpfile << " ECP gen " << endl;
#else
  inpfile << " BASIS 6-31G " << endl;     
#endif
//  inpfile << " WAVEFUNCTION_ANALYSIS FALSE " << endl;
#if DO_NBO
  inpfile << " NBO 1 " << endl;
  inpfile << " print_orbitals 250 " << endl;
  inpfile << " MOLDEN_FORMAT TRUE " << endl;
#endif
  inpfile << " $end " << endl;
  inpfile << endl;

#if (BASISMIX || BASISGEN)
  inpfile << endl;
  string cmd = "cat qmix >> "+inpfile_string;
  system(cmd.c_str());
  inpfile << endl;
#endif

  inpfile.close();

  return;

}

// for para execution
void DFT::ts_dnr(string filename) {

  //printf(" in DFT/ts_dnr() \n");

  energy0 = energy = 0;

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  inpfile << " $molecule " << endl;
  inpfile << " " << DFT_CHARGE << " " << DFT_SPIN << endl;

  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  inpfile << " $end " << endl;
  inpfile << endl << endl;
  inpfile << " $rem " << endl;
  inpfile << " JOBTYPE TS " << endl;
  inpfile << " GEOM_OPT_DMAX 80 " << endl;
  inpfile << " GEOM_OPT_MAX_CYCLES 150 " << endl;

  inpfile << " GEOM_OPT_TOL_DISPLACEMENT 2500 " << endl;
  inpfile << " GEOM_OPT_TOL_GRADIENT     800 " << endl;
  inpfile << " GEOM_OPT_TOL_ENERGY      5000 " << endl;

  inpfile << " EXCHANGE " << DFT_EXCHANGE << endl;
  inpfile << " CORRELATION " << DFT_CORRELATION << endl;
#if UNRESTRICTED
  inpfile << " UNRESTRICTED TRUE " << endl;
#endif
  inpfile << " SCF_ALGORITHM rca_diis " << endl;
//  inpfile << " SCF_ALGORITHM gdm " << endl;
  inpfile << " SCF_MAX_CYCLES 250 " << endl;
  inpfile << " SCF_CONVERGENCE 6 " << endl;
#if B631GSS
  inpfile << " BASIS 6-31G** " << endl;
#elif B631GS 
  inpfile << " BASIS 6-31G* " << endl;
#elif LANL2DZ
  inpfile << " BASIS LANL2DZ " << endl;
  inpfile << " ECP LANL2DZ " << endl;
#elif BASISMIX
  inpfile << " BASIS mixed " << endl;
  inpfile << " ECP gen " << endl;
#else
  inpfile << " BASIS 6-31G " << endl;     
#endif
//  inpfile << " WAVEFUNCTION_ANALYSIS FALSE " << endl;
#if DO_NBO
  inpfile << " NBO 1 " << endl;
  inpfile << " print_orbitals 250 " << endl;
#endif
  inpfile << " MOLDEN_FORMAT TRUE " << endl;
  inpfile << " $end " << endl;
  inpfile << endl;

#if (BASISMIX || BASISGEN)
  inpfile << endl;
  string cmd = "cat qmix >> "+inpfile_string;
  system(cmd.c_str());
  inpfile << endl;
#endif

  return;
}

void DFT::get_structure(string filename, double* xyzc) {

  xyz_read(filename);
  for (int i=0;i<natoms*3;i++)
    xyzc[i] = xyz[i];

  return;
}

double DFT::get_energy(string filename) {

  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening dft file: %s \n",oname.c_str()); return 0.0; }
  string line;
  vector<string> tok_line;
  energy = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (line.find("Total energy in the final basis set")!=string::npos)
    {
      cout << "  DFT out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[8].c_str());
      break;
    }
  }
 
//  printf(" DFT energy: %1.4f \n",energy); 

  if (abs(energy)<0.00001 || (energy != energy))
  {
    printf(" energy zero, DFT failed \n");
    return 10000;
  }

  return energy;
}


double DFT::get_opt_energy(string filename) {

  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening dft file: %s \n",oname.c_str()); return 0.0; }
  string line;
  vector<string> tok_line;
  energy = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (line.find("Total energy in the final basis set")!=string::npos)
    {
//      cout << "  DFT out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[8].c_str());
    }
  }
 
//  printf(" DFT opt energy: %1.4f \n",energy); 

  if (abs(energy)<0.00001 || (energy != energy))
  {
    printf(" energy zero, DFT failed \n");
    return 10000;
  }

  return energy;
}

double DFT::get_energy_ts(string filename) {

  energyts = 0.0;
  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening ts file: %s \n",oname.c_str()); return 0.0; }
  string line;
  vector<string> tok_line;
  int complete = 0;
  int hline = 0;
  double hval1 = -1;
  double gradts = -1;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;

    if (line.find("Gradient   ")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      gradts = atof(tok_line[1].c_str());
    }
    else if (line.find("Hessian Eigenvalues")!=string::npos)
    {
      hline++;
    }
    else if (hline)
    {
      tok_line = StringTools::tokenize(line, " \t");
      hval1 = atof(tok_line[0].c_str());
      hline--;
    }

    else if (line.find("Total energy in the final basis")!=string::npos)
    {
//      cout << "  DFT TS out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energyts=atof(tok_line[8].c_str());
    }
    else if (line.find("Final energy")!=string::npos)
    {
//      cout << "  DFT TS Final out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energyts=atof(tok_line[3].c_str());
      complete = 1;
      break;
    }
  }
 
//  printf(" DFT TS energy: %1.4f \n",energyts); 

  if (abs(energyts)<0.00001 || (energyts != energyts))
  {
    printf(" energy zero, DFT failed \n");
    return 10000;
  }
  if (!complete)
  {
    printf(" WARNING: %s not converged, energy is: %1.4f grad is: %1.4f eigenval: %1.3f \n",filename.c_str(),energyts,gradts,hval1);
    converged = 0;
  }
  else
  {
    printf(" SUCCESS: %s converged,     energy is: %1.4f grad is: %1.4f eigenval: %1.3f \n",filename.c_str(),energyts,gradts,hval1);
    converged = 1;
  }

  return energyts;
}

void DFT::alloc(int natoms_i) {
 
//  printf(" in DFT/alloc() \n");

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz = new double[3*natoms];

  return;
}

void DFT::init(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in DFT/init() \n");

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz = new double[3*natoms];

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz[i] = xyz_i[i];  

  return;
}

void DFT::freemem() {

  natoms = 0;
  delete [] xyz0;
  delete [] xyz;
  delete [] anumbers;
  delete [] anames;

  return;
}

void DFT::reset(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in DFT/reset() \n");
  natoms = natoms_i;

#if 0
  if (natoms!=natoms_i)
  {
    printf(" DFT reset failed due to different # of atoms \n");
    return;
  }
#endif

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz[i] = xyz_i[i];  

  return;
}

void DFT::xyz_read(string filename){ 
   
//  printf(" in DFT - xyz_read \n");
  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening dft xyz file: %s \n",oname.c_str()); return; }
  string line;
  vector<string> tok_line;
  int count = 0;

  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (count == 1)
    {
//      cout << " reading! " << endl;
      for (int i=0;i<natoms;i++)
      {
        getline(output,line);
        if (!StringTools::cleanstring(line))
        {
          count = 0;
          break;
        }
        if (line.find("----")!=string::npos)
        {
          printf(" ERROR: probably wrong # of atoms in dft output \n");
          printf(" stopped at atom %i \n",i+1);
          exit(1);
        }
        vector<string> tok_line = StringTools::tokenize(line, " \t");
        if (tok_line.size()<5)
        {
          printf("\n ERROR: cannot read DFT output \n");
          printf(" line: %s \n",line.c_str());
          return;
        }
        if (tok_line[1]!=anames[i])
        {
          printf("\n ERROR: anames does not match file output: %s vs %s \n",anames[i].c_str(),tok_line[1].c_str());
        }
        xyz[3*i+0]=atof(tok_line[2].c_str());
        xyz[3*i+1]=atof(tok_line[3].c_str());
        xyz[3*i+2]=atof(tok_line[4].c_str());
//      cout << tok_line[0] << " " << tok_line[1] << " " << tok_line[2] << " " << endl; 
      }
      count = 0;
    }

    if (line.find("Standard Nuclear Orientation")!=string::npos)
    {
//      cout << "Final: " << line << endl;
      count++;
      getline(output,line);
    }
  }

  output.close();
 
  return;
}   



void DFT::get_charges(string filename, double* q)
{ 
  for (int i=0;i<natoms;i++)
    q[i] = 0.;

//  printf(" in DFT - xyz_read \n");
  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening dft out file: %s \n",oname.c_str()); return; }
  string line;
  vector<string> tok_line;
  int count = 0;

  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (count == 1)
    {
//      cout << " reading! " << endl;
      for (int i=0;i<natoms;i++)
      {
        getline(output,line);
        if (!StringTools::cleanstring(line))
        {
          count = 0;
          break;
        }
        if (line.find("----")!=string::npos)
        {
          printf(" ERROR: probably wrong # of atoms in dft output \n");
          printf(" stopped at atom %i \n",i+1);
          exit(1);
        }
        vector<string> tok_line = StringTools::tokenize(line, " \t");
        if (tok_line.size()<3)
        {
          printf("\n ERROR: cannot read DFT output \n");
          printf(" line: %s \n",line.c_str());
          return;
        }
        q[i] = atof(tok_line[2].c_str());
      }
      count = 0;
    }

    if (line.find("Ground-State Mulliken Net Atomic Charges")!=string::npos)
    {
      count++;
      getline(output,line);
      getline(output,line);
    }
  }

  output.close();

#if 1
  printf("  charges:");
  for (int i=0;i<natoms;i++)
    printf(" %4.3f",q[i]);
  printf("\n");
#endif
#if 0
  for (int i=0;i<natoms;i++)
    q[i] = 0.;
  printf(" DEBUG: charges overwritten with zero \n");
#endif
 
  return;
}   


