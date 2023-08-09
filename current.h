#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
 

void get_current(Parameters &parameters, MatrixVectorType &gf_int_r, MatrixVectorType &gf_int_l,
    MatrixVectorType &gamma, double &current, int left_right);

void get_current_transmission(Parameters &parameters, MatrixVectorType &gf_int_r, 
    MatrixVectorType &gamma_l, MatrixVectorType &gamma_r, double &current);