#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include "aim.h"
#include <limits>
#include </usr/include/eigen3/Eigen/Dense>
#include "utilis.h"

AIM::AIM(const Parameters &parameters, const MatrixVectorType &local_gf_retarded, const MatrixVectorType &local_gf_lesser, int i, int spin)
{
    this->impurity_gf_mb_retarded.resize(parameters.steps_myid);
    this->dynamical_field_retarded.resize(parameters.steps_myid);
    this->self_energy_mb_retarded.resize(parameters.steps_myid);
    
    this->self_energy_mb_lesser.resize(parameters.steps_myid);
    this->impurity_gf_mb_lesser.resize(parameters.steps_myid);
    this->dynamical_field_lesser.resize(parameters.steps_myid);
    this->hybridisation_lesser.resize(parameters.steps_myid);
    
    this->fermi_function_eff.resize(parameters.steps_myid);

    get_dynamical_fields(parameters, local_gf_retarded, local_gf_lesser, i);
    //print_to_file(parameters, "dynamical_field_retarded", this->dynamical_field_retarded, 1);
    //print_to_file(parameters, "dynamical_field_lesser", this->dynamical_field_lesser, 1);

    if (parameters.impurity_solver == 1) {//brute force method
        if (parameters.read_se == 1) {
            std::cout << "Reading the many-body self energy from file \n";
            read_self_energy(parameters, spin);
        }

        get_impurity_gf_mb_retarded(parameters);
        get_lesser_hybridisation(parameters);   
        get_impurity_gf_mb_lesser(parameters);
        print_to_file(parameters, "hybridisation_lesser", this->hybridisation_lesser, 1);

    } else if (parameters.impurity_solver == 2){
        get_effective_fermi_function(parameters); 
        //print_to_file(parameters, "effective_fermi_function", this->fermi_function_eff, 1);
    }
    
    
    //for (int r = 0; r < parameters.steps_myid; r++) std::cout << parameters.energy[r] << " " << local_gf_retarded[r](0, 0) << " " << this->dynamical_field_retarded[r] 
    //    << " " << this->impurity_gf_mb_retarded[r] << "\n";
}

void AIM::get_dynamical_fields(const Parameters &parameters, const MatrixVectorType &local_gf_retarded, 
    const MatrixVectorType &local_gf_lesser, int i)
{
    for (int r = 0; r < parameters.steps_myid; r++){
        this->dynamical_field_retarded[r] = local_gf_retarded[r](i, i);
        this->dynamical_field_lesser[r] = local_gf_lesser[r](i, i);
    }
}

void AIM::get_lesser_hybridisation(const Parameters &parameters)
{  
    for (int r = 0; r < parameters.steps_myid; r++){
        this->hybridisation_lesser[r] = (1.0 / this->dynamical_field_retarded[r]) * this->dynamical_field_lesser[r] *
                                           1.0 / std::conj(this->dynamical_field_retarded[r]);
    }
}     

void AIM::get_impurity_gf_mb_retarded(const Parameters &parameters)
{
    for (int r = 0; r < parameters.steps_myid; r++){
        this->impurity_gf_mb_retarded[r] = 1.0 / ((1.0 / this->dynamical_field_retarded[r]) - this->self_energy_mb_retarded[r]);
        //std::cout << this->impurity_gf_mb_retarded[r] << " " << this->dynamical_field_retarded[r] << " \n";
    }
    //std::cout << std::setprecision(15) << this->impurity_gf_mb_retarded[0] << " " << this->self_energy_mb_retarded[0] << " " << this->dynamical_field_retarded[0] << " \n";

}

void AIM::get_impurity_gf_mb_lesser(const Parameters &parameters)
{
    for (int r = 0; r < parameters.steps_myid; r++){
        dcomp advanced = std::conj(this->impurity_gf_mb_retarded[r]);
        this->impurity_gf_mb_lesser[r] = (this->impurity_gf_mb_retarded[r] * (this->hybridisation_lesser[r]
            + this->self_energy_mb_lesser[r]) * advanced); //need to be sure that this is indeed imaginary
        
        //std::cout << impurity_gf_mb_lesser[r] <<  "  " << this->impurity_gf_mb_retarded[r] << " " << " " << 
        //    this->hybridisation_lesser[r] << " " << this->self_energy_mb_lesser[r] << " " << advanced  << "  " << r << std::endl;
    }
}

void AIM::get_effective_fermi_function(const Parameters &parameters){
    //MatrixVectorType fermi_function(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));
    for (int r = 0; r < parameters.steps_myid; r++){
        dcomp advanced = std::conj(this->dynamical_field_retarded[r]);
        this->fermi_function_eff[r] =  - this->dynamical_field_lesser[r].imag() / ( (this->dynamical_field_retarded[r] - advanced).imag());
        //std::cout << this->fermi_function_eff[r] << " " << this->dynamical_field_lesser[r].imag() << " " << (this->dynamical_field_retarded[r] - advanced) << std::endl;
        //fermi_function[r](0, 0) = this->fermi_function_eff[r];
    }
    //print_to_file(parameters, "fermi_function", fermi_function, 0, 0, 1);
}

void AIM::read_self_energy(const Parameters &parameters, int spin) {
    std::string var_retarded, var_lesser;
    if (spin == 1) {// 1 is for spin up
        var_lesser = "1_se_lesser_up_0_0.dat";
        var_retarded = "1_se_retarded_up_0_0.dat";
    } else if (spin == 2) {//2 is for spin down
        var_lesser = "1_se_lesser_down_0_0.dat";
        var_retarded = "1_se_retarded_down_0_0.dat";
    }

    fstream my_file_retarded;
	my_file_retarded.open(var_retarded, ios::in);
	std::string line;
    int lineCount = -1;
    int energy_index = 0;
    while (std::getline(my_file_retarded, line)) {
        lineCount++;
        if ((energy_index >= parameters.start[parameters.myid]) && (energy_index < parameters.start[parameters.myid] + parameters.steps_myid)) {
            //this selects the appriotate lines for each process
            std::istringstream iss(line);
            double num1, num2, num3;
            if (iss >> num1 >> num2 >> num3) {
                // Use the values stored in num1, num2, and num3
                int y = energy_index - parameters.start[parameters.myid];
                this->self_energy_mb_retarded[y] = num2 + parameters.j1 * num3;
            } else {
                std::cout << "Invalid line format: " << line << std::endl;
            }
        }
        energy_index++;        
    }
    my_file_retarded.close();

    fstream my_file_lesser;
	my_file_lesser.open(var_lesser, ios::in);
    lineCount = -1;
    energy_index = 0;
    while (std::getline(my_file_lesser, line)) {
        lineCount++;
        if ((energy_index >= parameters.start[parameters.myid]) && (energy_index < parameters.start[parameters.myid] + parameters.steps_myid)) {
            //this selects the appriotate lines for each process
            std::istringstream iss(line);
            double num1, num2, num3;
            if (iss >> num1 >> num2 >> num3) {
                // Use the values stored in num1, num2, and num3
                int y = energy_index - parameters.start[parameters.myid];
                this->self_energy_mb_lesser[y] = num2 + parameters.j1 * num3;
            } else {
                std::cout << "Invalid line format: " << line << std::endl;
            }
        }
        energy_index++;        
    }
    my_file_lesser.close();

}
