#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
#include "current.h"
#include <mpi.h>
using namespace std;

void get_current(Parameters &parameters, MatrixVectorType &gf_int_r, MatrixVectorType &gf_int_l,
    MatrixVectorType &gamma, double &current, int left_right) {
        dcomp integrand_myid = 0.0, integrand = 0.0;
        if (left_right == 0) {
            for (int r = 0; r < parameters.steps_myid; r++) {
                int y = r + parameters.start.at(parameters.myid);
                Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
                Eigen::MatrixXcd spectral_function = parameters.j1 * (gf_int_r.at(r) - gf_int_a);
                Eigen::MatrixXcd trace =  fermi_function(parameters.energy.at(y) - parameters.voltage, parameters) *
                        gamma.at(r) * spectral_function + parameters.j1 * gamma.at(r) * gf_int_l.at(r);

                integrand_myid += trace.trace();
            }
            MPI_Allreduce(&integrand_myid, &integrand, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        } else {
            for (int r = 0; r < parameters.steps_myid; r++) {
                int y = r + parameters.start.at(parameters.myid);
                Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
                Eigen::MatrixXcd spectral_function = parameters.j1 * (gf_int_r.at(r) - gf_int_a);
                Eigen::MatrixXcd trace =  fermi_function(parameters.energy.at(y) + parameters.voltage, parameters) *
                        gamma.at(r) * spectral_function + parameters.j1 * gamma.at(r) * gf_int_l.at(r);
                integrand_myid += trace.trace();
            }
            MPI_Allreduce(&integrand_myid, &integrand, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        }
        //cout << integrand << endl;
        current = (integrand.real()) * (parameters.energy.at(1) - parameters.energy.at(0)) * 0.00005270602546 * 0.5;
}


//void get_current_transmission(Parameters &parameters, MatrixVectorType &gf_int_r, 
//    MatrixVectorType &gamma_l, MatrixVectorType &gamma_r, double &current) {
//    double integrand_myid = 0, integrand = 0;
//
//    std::ofstream my_file;
//    std::ostringstream oss_right;
//    oss_right << parameters.path << "/transmission.dat";
//    std::string var = oss_right.str();
//    //std::cout << var_right << std::endl;
//	my_file.open(var);
//    for (int r =0; r < parameters.steps_myid; r++) {
//        int y = r + parameters.start.at(parameters.myid);
//        Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint();
//        dcomp transmission = (gamma_l.at(r) * gf_int_r.at(r) * gamma_r.at(r) * gf_int_a).trace();
//        my_file << parameters.energy.at(r) << "  " << transmission.real() << "  " << transmission.imag() << "\n";
//        integrand_myid += transmission.real() * (fermi_function(parameters.energy.at(y) - parameters.voltage, parameters) - 
//            fermi_function(parameters.energy.at(y) + parameters.voltage, parameters));
//    }
//    my_file.close();
//    MPI_Allreduce(&integrand_myid, &integrand, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
//    current = 0.5 * (integrand) * (parameters.energy.at(1) - parameters.energy.at(0)) * 0.00005270602546;
//}


void get_current_transmission(Parameters &parameters, MatrixVectorType &gf_int_r, 
    MatrixVectorType &gamma_l, MatrixVectorType &gamma_r, double &current) {
    double integrand_myid = 0, integrand = 0;

    std::ofstream my_file;
    std::ostringstream oss_right;
    oss_right << parameters.path << "/transmission.dat";
    std::string var = oss_right.str();
    //std::cout << var_right << std::endl;
    my_file.open(var);
    for (int r =0; r < parameters.steps_myid; r++) {
        int y = r + parameters.start.at(parameters.myid);
        Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint();
        dcomp transmission = (gamma_l.at(r) * gf_int_r.at(r) * gamma_r.at(r) * gf_int_a).trace();
        my_file << parameters.energy.at(r) << "  " << transmission.real() << "  " << transmission.imag() << "\n";
        integrand_myid += transmission.real() * (fermi_function(parameters.energy.at(y) - parameters.voltage, parameters) -
            fermi_function(parameters.energy.at(y) + parameters.voltage, parameters));
    }
    my_file.close();
    MPI_Allreduce(&integrand_myid, &integrand, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    current = 0.5 * (integrand) * (parameters.energy.at(1) - parameters.energy.at(0)) * 0.00005270602546;
}

