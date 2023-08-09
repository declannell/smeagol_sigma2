#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>

void get_number_energy_points(Parameters &parameters);

void read_non_interacting_gf(Parameters &parameters, MatrixVectorType &gf_non_int_r, int spin);

void read_gamma(Parameters &parameters, MatrixVectorType &gamma_left, MatrixVectorType &gamma_right, int spin);

void get_lesser_se(Parameters &parameters, MatrixVectorType &se_lesser, MatrixVectorType &sgamma, int left_right);

void get_lesser_gf(Parameters &parameters, MatrixVectorType &gf_int_r, MatrixVectorType &se_left,
     MatrixVectorType &se_right, MatrixVectorType &sigma_mb_l, MatrixVectorType &gf_int_l);

void get_interacting_retarded_gf(const Parameters &parameters, MatrixVectorType &gf_int_r, MatrixType &hamiltonian,
    const MatrixVectorType &sigma_mb_r, const MatrixVectorType &embedding_r, const MatrixVectorType &embedding_l);

void read_hamiltonian(const Parameters &parameters, MatrixType &hamiltonian_up, MatrixType &hamiltonian_down);

void get_retarded_se(Parameters &parameters, MatrixVectorType &se_retarded, MatrixVectorType &gamma);

void read_delta(Parameters &parameters, MatrixVectorType &delta, int spin);

void get_interacting_retarded_gf(const Parameters &parameters, MatrixVectorType &gf_int_r, MatrixType &hamiltonian,
    const MatrixVectorType &sigma_mb_r, const MatrixVectorType &delta);
    
void get_model_gf(const Parameters &parameters, MatrixVectorType &gf_int_r, const MatrixVectorType &sigma_mb_r,
    const MatrixVectorType &embedding_r, const MatrixVectorType &embedding_l);
    

