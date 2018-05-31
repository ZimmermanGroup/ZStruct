#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>

#include <sys/time.h>
#include "stringtools.h"

double getRa(int* anumbers, int a1);

void sort_2(double* v1);
void sort_3(double* v1);
void sort_4(double* v1);
void sort_5(double* v1);
void sort_2p(double* v1, double* v2);
void sort_3p(double* v1, double* v2);
void sort_4p(double* v1, double* v2);
void sort_5p(double* v1, double* v2);
void sort_n2p(double* v1, double* v2);
void sort_n3p(double* v1, double* v2);
void sort_n4p(double* v1, double* v2);
void sort_n5p(double* v1, double* v2);

int sign(double x);
void cross(double* x, double* r1, double* r2);
double norm(double* x, int size);
int randomi(int a);
double randomf(double a);
int close_val(double a, double b, double c);
int fact(int n);
double sigmoid(double A, double A0, double alpha);

void trans(double* A, double* B, int m, int n);
int Diagonalize(double* A, double* eigen, int size);
int Invert(double* A, int m);
int matmult_sq(double* A, double* B, double* C, int N);

void order_array(int* olist, double* Ao, double* A, int na);
void get_rotation_matrix(double** rotMat, double* thetas);

int check_array(int size, double* A);
double read_xyz_file(string filename, int natoms, double* xyz1);
double read_chk_file(string filename, int natoms, double* xyz1);
