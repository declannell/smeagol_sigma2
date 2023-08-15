#pragma once
#include <complex>
#include <cmath>
#include <vector>
#include </usr/include/eigen3/Eigen/Dense>

using MatrixType = Eigen::MatrixXcd;
using MatrixVectorType = std::vector<MatrixType>;


typedef std::complex<double> dcomp;


struct Parameters
{
    double chemical_potential;
    int steps; // number of energy points we take
    std::vector<double> energy;
    dcomp j1; // this is a complex number class defined within the complex library
    int num_orbitals;
    double voltage;
    static Parameters from_file();
    int interaction;
    int system; //if system ==0, the code is for Fe_MgO. If system == 1 the code is for the copper system. 
    std::string path;
    int impurity_solver; //1 is brute force. 2 is kramer-kronig
    int size; //the size of the communicator
    int myid; //the id of the process
    std::vector<int> start; //the starting index of the energy array for each process
    std::vector<int> end; //the ending index of the energy array for each process
    int steps_myid; //this is the number of steps the process has
    std::vector<int> steps_proc; //this is the number of steps the other processes have
    std::vector<int> displs;
    double delta_energy;
    double hubbard_interaction;
    int self_consistent_steps;
    double convergence;
    int grid_density;
    int read_gf; //1 means non-interacting gf read from file. 0 means it is calculated via G = (E- H -sigma)^{-1}
    double gamma; //broadening of the leads
    double onsite;
    int model_calc;//1 means we do an ab initio calculation. 0 means we use the model parameters.
    double random_err;
    int num_atoms; //total number of atoms in the scattering region. 
    int num_orb_total; //total number of orbtial in the scattering region
    double se_mixing;// this is a number between 1 and zero which determines the mixing of the self energy between each iteration
};




double fermi_function(double energy, const Parameters &parameters);
void print_parameters(Parameters& parameters);
