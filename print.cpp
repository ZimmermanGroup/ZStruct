// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //


#include "icoord.h"

void ICoord::print_xyz(){

 printf(" %i \n",natoms);
 printf("\n");
 for (int i=0;i<natoms;i++) 
 {
     cout << "  " << anames[i];
     printf(" %f %f %f \n",coords[3*i+0],coords[3*i+1],coords[3*i+2]);
 }
// printf("\n");

}

void ICoord::print_ic()
{
   printf("\n printing internals \n");
   printf(" number of bonds: %i\n",nbonds);
   for (int i=0;i<nbonds;i++)
      printf(" bond %2i: %2i to %2i: %1.3f \n",i+1,bonds[i][0]+1,bonds[i][1]+1,bondd[i]);
   printf("\n");
  
   for (int i=0;i<natoms;i++)
      printf(" atom %2i is %i coordinate (%s) \n",i+1,coordn[i],anames[i].c_str());
   printf("\n");


   printf(" number of angles: %2i\n",nangles);
   for (int i=0;i<nangles;i++)
      printf(" angle %2i: %2i %2i %2i: %1.1f \n",i+1,angles[i][0]+1,angles[i][1]+1,angles[i][2]+1,anglev[i]);
   printf("\n");

#if 0
   printf(" number of torsions: %2i \n",ntor);
   for (int i=0;i<ntor;i++)
      printf(" torsion %2i: %2i %2i %2i %2i: %1.1f \n",i+1,torsions[i][0]+1,torsions[i][1]+1,torsions[i][2]+1,torsions[i][3]+1,torv[i]);
   printf("\n");
#endif

#if 1
   printf(" number of improper torsions: %2i \n",nimptor);
   for (int i=0;i<nimptor;i++)
      printf(" imptor %i: %2i %2i %2i %2i: %1.1f \n",i+1,imptor[i][0]+1,imptor[i][1]+1,imptor[i][2]+1,imptor[i][3]+1,imptorv[i]);
   printf("\n");
#endif

   printf(" number of nonbonds: %3i \n",n_nonbond);
#if 0
   for (int i=0;i<n_nonbond;i++)
      printf(" nonbond %2i: %2i to %2i: %1.2f \n",i+1,nonbond[i][0]+1,nonbond[i][1]+1,nonbondd[i]);
   printf("\n");
#endif
   printf("\n");
}

void ICoord::print_bonds()
{
   printf("\n printing internals \n");
   printf(" number of bonds: %2i \n",nbonds);
   for (int i=0;i<nbonds;i++)
      printf(" bond %2i: %2i to %2i: %1.4f \n",i+1,bonds[i][0]+1,bonds[i][1]+1,bondd[i]);
   printf("\n");
  
   for (int i=0;i<natoms;i++)
      printf(" atom %2i is %2i coordinate (%s) \n",i+1,coordn[i],anames[i].c_str());
   printf("\n");
}


void print_xyz_gen(int natoms, string* anames, double* coords){
   
 printf(" %i \n",natoms);
 printf("\n");
 for (int i=0;i<natoms;i++)
 {
     cout << "  " << anames[i];
     printf(" %f %f %f \n",coords[3*i+0],coords[3*i+1],coords[3*i+2]);
 }
// printf("\n");
   
}

void print_triple_xyz(int natoms, string* anames, int* anumbers, double* coords0, double* coords1, double* coords2, double* e){
   
 printf(" %i \n",natoms);
 printf("\n");
 for (int i=0;i<natoms;i++)
 {
     cout << "  " << anames[i];
     printf(" %f %f %f \n",coords0[3*i+0],coords0[3*i+1],coords0[3*i+2]);
 }
 printf(" %i \n",natoms);
 printf("\n");
 for (int i=0;i<natoms;i++)
 {
     cout << "  " << anames[i];
     printf(" %f %f %f \n",coords1[3*i+0],coords1[3*i+1],coords1[3*i+2]);
 }
 printf(" %i \n",natoms);
 printf("\n");
 for (int i=0;i<natoms;i++)
 {
     cout << "  " << anames[i];
     printf(" %f %f %f \n",coords2[3*i+0],coords2[3*i+1],coords2[3*i+2]);
 }
// printf("\n");
   
}

void ICoord::print_xyz_save(string xyzfile_string){

  ofstream xyzfile;
//  string xyzfile_string = "xyzfile.txt";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(8);

  int natoms1 = natoms;
  for (int i=0;i<natoms1;i++)
  if (anames[i]=="X")
    natoms1--;

  xyzfile << " " << natoms1 << endl;
  xyzfile << "  " << q1 << endl;
  for (int i=0;i<natoms;i++) 
  if (anames[i]!="X")
  {
    xyzfile << "  " << anames[i];
    xyzfile << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2];
    xyzfile << endl;
  }
  //xyzfile << endl;
  
  xyzfile.close();

  return;
}

void ICoord::print_xyz_save_no_natoms(string xyzfile_string){

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(8);

  int natoms1 = natoms;
  for (int i=0;i<natoms1;i++)
  if (anames[i]=="X")
    natoms1--;

  for (int i=0;i<natoms;i++) 
  if (anames[i]!="X")
  {
    xyzfile << "  " << anames[i];
    xyzfile << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2];
    xyzfile << endl;
  }
  
  xyzfile.close();

  return;
}

void ICoord::print_xyz_save(string xyzfile_string, double energy){

  ofstream xyzfile;
//  string xyzfile_string = "xyzfile.txt";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(6);

  int natoms1 = natoms;
  for (int i=0;i<natoms1;i++)
  if (anames[i]=="X")
    natoms1--;

   xyzfile << " " << natoms1 << endl;
   xyzfile << " " << energy << endl;
   for (int i=0;i<natoms;i++) 
   if (anames[i]!="X")
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2];
     xyzfile << endl;
   }

  xyzfile.close();

  return;
}


void print_single_xyz_save(string xyzfile_string, int natoms, string* anames, double* coords0){

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(6);

  int natoms1 = natoms;
  for (int i=0;i<natoms1;i++)
  if (anames[i]=="X")
    natoms1--;

   xyzfile << " " << natoms1 << endl;
   xyzfile << endl;
   for (int i=0;i<natoms;i++) 
   if (anames[i]!="X")
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << coords0[3*i+0] << " " << coords0[3*i+1] << " " << coords0[3*i+2];
     xyzfile << endl;
   }

  xyzfile.close();

  return;
}

void print_double_xyz_save(string xyzfile_string, int natoms, string* anames, double* coords0, double* coords1, double* e){

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(6);

  int natoms1 = natoms;
  for (int i=0;i<natoms1;i++)
  if (anames[i]=="X")
    natoms1--;

   xyzfile << " " << natoms1 << endl;
   xyzfile << " " << e[0] << endl;
   for (int i=0;i<natoms;i++) 
   if (anames[i]!="X")
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << coords0[3*i+0] << " " << coords0[3*i+1] << " " << coords0[3*i+2];
     xyzfile << endl;
   }
   xyzfile << " " << natoms1 << endl;
   xyzfile << " " << e[1] << endl;
   for (int i=0;i<natoms;i++) 
   if (anames[i]!="X")
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << coords1[3*i+0] << " " << coords1[3*i+1] << " " << coords1[3*i+2];
     xyzfile << endl;
   }
   xyzfile << endl;

  xyzfile.close();

  return;
}

void print_triple_xyz_save(string xyzfile_string, int natoms, string* anames, double* coords0, double* coords1, double* coords2, double* e){

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile << setprecision(6);

  int natoms1 = natoms;
  for (int i=0;i<natoms1;i++)
  if (anames[i]=="X")
    natoms1--;

   xyzfile << " " << natoms1 << endl;
   xyzfile << " " << e[0] << endl;
   for (int i=0;i<natoms;i++) 
   if (anames[i]!="X")
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << coords0[3*i+0] << " " << coords0[3*i+1] << " " << coords0[3*i+2];
     xyzfile << endl;
   }
   xyzfile << " " << natoms1 << endl;
   xyzfile << " " << e[1] << endl;
   for (int i=0;i<natoms;i++) 
   if (anames[i]!="X")
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << coords1[3*i+0] << " " << coords1[3*i+1] << " " << coords1[3*i+2];
     xyzfile << endl;
   }
   xyzfile << " " << natoms1 << endl;
   xyzfile << " " << e[2] << endl;
   for (int i=0;i<natoms;i++) 
   if (anames[i]!="X")
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << coords2[3*i+0] << " " << coords2[3*i+1] << " " << coords2[3*i+2];
     xyzfile << endl;
   }

  xyzfile.close();

  return;
}


void save_gephi(int npts1, int* ids, double* distances, double* values, double alpha, double savethresh)
{
  ofstream gfile;
  string gfile_string = "1kgraph_n.csv";
  gfile.open(gfile_string.c_str());
  gfile.setf(ios::fixed);
  gfile.setf(ios::left);
  gfile << setprecision(3);

  gfile << "Id;Label;Pr" << endl;
  for (int i=0;i<npts1;i++)
  if (values[i]!=111.1 && ids[i]<999999)
    gfile << i << ";" << ids[i] << ";" << values[i] << endl;
  gfile.close();
    
  gfile_string = "1kgraph_e.csv";
  gfile.open(gfile_string.c_str());
  gfile.setf(ios::fixed);
  gfile.setf(ios::left);
  gfile << setprecision(3);

  gfile << "Source;Target;Id;Weight" << endl;
  for (int i=0;i<npts1;i++)
  if (values[i]!=111.1 && ids[i]<999999)
  {
    for (int j=0;j<npts1;j++)
    if (values[j]!=111.1 && ids[j]<999999)
    if (i!=j)
    {
      double val = exp(-distances[i*npts1+j]/alpha);
      //double val = distances[i*npts1+j];
      if (val>savethresh)
        gfile << i << ";" << j << ";" << i << ";" << val << endl;
    }
  }
  gfile.close();

#if 0
  printf("Id;label;Pr\n");
  for (int n=0;n<npts1;n++)
    printf(" %i;%i;%5.4f\n",n,n,values[n]);

  printf("Source;Target;Id;Weight\n");
  for (int n=0;n<npts1;n++)
  {
    for (int m=0;m<n;m++)
      printf(" %i;%i;%i;%5.4f\n",n,m,n,exp(-distances[n*npts1+m]/KNNR_ALPHA));
  }
#endif

  return;
} 

void save_gephi_2(int npts1, int* ids, double* distances, double* values, double* values_print, double alpha, double savethresh)
{
  ofstream gfile;
  string gfile_string = "1kgraph_n.csv";
  gfile.open(gfile_string.c_str());
  gfile.setf(ios::fixed);
  gfile.setf(ios::left);
  gfile << setprecision(3);

  gfile << "Id;Label;Pr;PrT2" << endl;
  for (int i=0;i<npts1;i++)
  if (values[i]!=111.1 && ids[i]<999999)
    gfile << i << ";" << ids[i] << ";" << values_print[i] << ";" << values[i] << endl;
  gfile.close();
    
  gfile_string = "1kgraph_e.csv";
  gfile.open(gfile_string.c_str());
  gfile.setf(ios::fixed);
  gfile.setf(ios::left);
  gfile << setprecision(3);

  gfile << "Source;Target;Id;Weight" << endl;
  for (int i=0;i<npts1;i++)
  if (values[i]!=111.1 && ids[i]<999999)
  {
    for (int j=0;j<npts1;j++)
    if (values[j]!=111.1 && ids[j]<999999)
    if (i!=j)
    {
      double val = exp(-distances[i*npts1+j]/alpha);
      //double val = distances[i*npts1+j];
      if (val>savethresh)
        gfile << i << ";" << j << ";" << i << ";" << val << endl;
    }
  }
  gfile.close();

#if 0
  printf("Id;label;Pr\n");
  for (int n=0;n<npts1;n++)
    printf(" %i;%i;%5.4f\n",n,n,values[n]);

  printf("Source;Target;Id;Weight\n");
  for (int n=0;n<npts1;n++)
  {
    for (int m=0;m<n;m++)
      printf(" %i;%i;%i;%5.4f\n",n,m,n,exp(-distances[n*npts1+m]/KNNR_ALPHA));
  }
#endif

  return;
} 

void save_gephi_3(int npts1, int* ids, double* distances, double* values, double* values_print, double alpha, double savethresh)
{
  printf("  in save_gephi_3, writing knn.csv files \n");
  int* saveme = new int[npts1*npts1];
  for (int i=0;i<npts1*npts1;i++)
    saveme[i] = 0;

  for (int i=0;i<npts1;i++)
  if (values[i]!=111.1 && ids[i]<999999)
  for (int j=0;j<npts1;j++)
  if (j!=i)
  if (values[j]!=111.1 && ids[j]<999999)
  {
    double val = exp(-distances[i*npts1+j]/alpha);
    if (val>savethresh)
      saveme[i*npts1+j] = saveme[j*npts1+i] = 1;
  }

  for (int i=0;i<npts1;i++)
  if (values[i]!=111.1 && ids[i]<999999)
  {
    int found = 0;
    for (int j=0;j<npts1;j++)
    if (saveme[i*npts1+j])
    { 
      found = 1;
      break;
    }

    if (!found)
    {
      double vmax = 0;
      int nmax = -1;
      for (int j=0;j<npts1;j++)
      if (j!=i)
      if (values[j]!=111.1 && ids[j]<999999)
      {
        double val = exp(-distances[i*npts1+j]/alpha);
        if (val>vmax)
        {
          nmax = j;
          vmax = val;
        }
      }
      if (nmax<0)
      {
        printf(" ERROR: nmax: %i  \n",nmax); 
        exit(1);
      }
      double val = exp(-distances[i*npts1+nmax]/alpha);
      if (val>savethresh/10.)
        saveme[i*npts1+nmax] = saveme[nmax*npts1+i] = 1;
    } 
  }

  ofstream gfile;
  string gfile_string = "1kgraph_n.csv";
  gfile.open(gfile_string.c_str());
  gfile.setf(ios::fixed);
  gfile.setf(ios::left);
  gfile << setprecision(3);

  gfile << "Id;Label;Pr;PrT2" << endl;
  for (int i=0;i<npts1;i++)
  if (values[i]!=111.1 && ids[i]<999999)
    gfile << i << ";" << ids[i] << ";" << values_print[i] << ";" << values[i] << endl;
  gfile.close();
    
  gfile_string = "1kgraph_e.csv";
  gfile.open(gfile_string.c_str());
  gfile.setf(ios::fixed);
  gfile.setf(ios::left);
  gfile << setprecision(3);

  gfile << "Source;Target;Id;Weight" << endl;
  for (int i=0;i<npts1;i++)
  {
    for (int j=0;j<i;j++)
    if (saveme[i*npts1+j])
    {
      double val = exp(-distances[i*npts1+j]/alpha);
      gfile << i << ";" << j << ";" << i << ";" << val << endl;
    }
  }
  gfile.close();

#if 0
  printf("Id;label;Pr\n");
  for (int n=0;n<npts1;n++)
    printf(" %i;%i;%5.4f\n",n,n,values[n]);

  printf("Source;Target;Id;Weight\n");
  for (int n=0;n<npts1;n++)
  {
    for (int m=0;m<n;m++)
      printf(" %i;%i;%i;%5.4f\n",n,m,n,exp(-distances[n*npts1+m]/KNNR_ALPHA));
  }
#endif

  delete [] saveme;

  return;
} 
