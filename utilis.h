#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip> 
#include "parameters.h"
using namespace std;

void decomp(int steps, int size, int myid, int *s, int *e);

void print_to_file(const Parameters &parameters, string filename, MatrixVectorType &quantity, int orbital_1, int orbital_2, int voltage_step);

double kramer_kronig_relation(const Parameters &parameters, std::vector<double> &impurity_self_energy_imag, int r);

double absolute_value(double num1);

void distribute_to_procs(const Parameters &parameters, std::vector<dcomp> &vec_1, std::vector<dcomp> &vec_2);

void distribute_to_procs(const Parameters &parameters, std::vector<double> &vec_1, const std::vector<double> &vec_2);

MatrixVectorType initializeMatrixVector(int size, int rows, int cols);

void print_to_file(const Parameters &parameters, string filename, std::vector<double> &quantity, int voltage_step);

void print_to_file(const Parameters &parameters, string filename, std::vector<dcomp> &quantity, int voltage_step);

