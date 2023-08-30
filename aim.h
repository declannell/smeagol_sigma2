#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include </usr/include/eigen3/Eigen/Dense>

typedef std::complex<double> dcomp;

class AIM
{
public:
    std::vector<dcomp> impurity_gf_mb_retarded;
    std::vector<dcomp> dynamical_field_retarded;
    std::vector<dcomp> self_energy_mb_retarded;
    std::vector<dcomp> self_energy_mb_lesser;
    std::vector<dcomp> impurity_gf_mb_lesser;
    std::vector<dcomp> hybridisation_lesser;
    std::vector<dcomp> dynamical_field_lesser;
    
    std::vector<double> fermi_function_eff;

    AIM(const Parameters &parameters, const MatrixVectorType &local_gf_retarded, const MatrixVectorType &local_gf_lesser, int i, int spin); 

    void get_dynamical_fields(const Parameters &parameters, const MatrixVectorType &local_gf_retarded,
         const MatrixVectorType &local_gf_lesser, int i);

    void get_impurity_gf_mb_retarded(const Parameters &parameters);

    void get_effective_fermi_function(const Parameters &parameters);

    void get_lesser_hybridisation(const Parameters &parameters);

    void get_impurity_gf_mb_lesser(const Parameters &parameters);
    void read_self_energy(const Parameters &parameters, int spin);
};
