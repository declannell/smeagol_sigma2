#include "parameters.h"
#include "gf.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
#include "utilis.h"
#include <mpi.h>
using namespace std;


void get_number_energy_points(Parameters &parameters) {
    fstream my_file;
    std::ostringstream oss;
    oss << parameters.path << "/Av-k_ReImGF_ReEne_K_1_1_1.dat";
    std::string var = oss.str();
    //std::cout << var << std::endl;
	my_file.open(var, ios::in);
	std::string line;
    parameters.steps = 0;
    //cout << gf_non_int_r.at(i).size() << endl;
    if (my_file.is_open()) {
        std::string line;
        int lineCount = -1;
        while (std::getline(my_file, line)) {
            lineCount++;
            if (lineCount % parameters.grid_density == 0) {//this selects every 10th lines
                std::istringstream iss(line);
                double num1, num2, num3;
                if (iss >> num1 >> num2 >> num3) {
                    // Use the values stored in num1, num2, and num3
                    parameters.energy.push_back(num1);
                } 
            }
        }
        parameters.steps = parameters.energy.size();
    } else {
        std::cout << var << " File cannot be opened for reading." << std::endl;
    }
    //std::cout << parameters.steps << std::endl;
    my_file.close(); 
}


void read_non_interacting_gf(Parameters &parameters, MatrixVectorType &gf_non_int_r, int spin) {
        for (int i = 0; i < parameters.num_orb_total; i++) {
            for (int j = 0; j < parameters.num_orb_total; j++) {
                fstream my_file;
                std::ostringstream oss;
                oss << parameters.path << "/Av-k_ReImGF_ReEne_K_" << i + 1 <<  "_" << j + 1 << "_" << spin << ".dat";
                std::string var = oss.str();
                //std::cout << var_left << std::endl;
	            my_file.open(var, ios::in);
	    	    std::string line;
                std::ifstream file(var); // Replace "input.txt" with the name of your input file

                int lineCount = -1;
                int energy_index = 0;
                while (std::getline(my_file, line)) {
                    lineCount++;
                    if (lineCount % parameters.grid_density == 0) {//this selects every 10th lines
                        if ((energy_index >= parameters.start.at(parameters.myid)) && (energy_index < parameters.start.at(parameters.myid) + parameters.steps_myid)) {
                            //this selects the appriotate lines for each process
                            std::istringstream iss(line);
                            double num1, num2, num3;
                            if (iss >> num1 >> num2 >> num3) {
                                // Use the values stored in num1, num2, and num3
                                int y = energy_index - parameters.start.at(parameters.myid);
                                gf_non_int_r.at(y)(i, j) = num2 + parameters.j1 * num3;

                            } else {
                                std::cout << "Invalid line format: " << line << std::endl;
                            }
                        }

                        energy_index++;
                    }
                }
                my_file.close();
            }
        }
}

void read_delta(Parameters &parameters, MatrixVectorType &delta, int spin) {

    fstream my_file;
    std::ostringstream oss;
    oss << parameters.path << "/Av-k_ReImDelta_ReEne_K_1_1_" << spin << ".dat";
    std::string var = oss.str();

    std::ifstream file(var); // Replace "input.txt" with the name of your input file
    if (!file.is_open()) {
        std::cout << "Failed to open the file. " << var << std::endl;
    }

    std::string line;
    int lineCount = -1;
    int energy_index = 0;
    while (std::getline(file, line)) {
        lineCount++;
        if (lineCount % parameters.grid_density == 0) {//this selects every 10th lines
            if ((energy_index >= parameters.start.at(parameters.myid)) && (energy_index < parameters.start.at(parameters.myid) + parameters.steps_myid)) {
                //this selects the appriotate lines for each process
                std::istringstream iss(line);
                double num1, num2, num3;
                if (iss >> num1 >> num2 >> num3) {
                    // Use the values stored in num1, num2, and num3
                    int y = energy_index - parameters.start.at(parameters.myid);
                    delta.at(y)(0, 0) = num2 + parameters.j1 * num3;
                } else {
                    std::cout << "Invalid line format: " << line << std::endl;
                }
            }

            energy_index++;
        }
    }

    file.close();
}

void read_gamma(Parameters &parameters, MatrixVectorType &gamma_left, MatrixVectorType &gamma_right, int spin) {
    if (parameters.model_calc == 1) { //this does the model ab initio calculation
        for (int i = 0; i < parameters.num_orb_total; i++) {
            for (int j = 0; j < parameters.num_orb_total; j++) {

                fstream my_file_left;
                std::ostringstream oss_left;
                oss_left << parameters.path << "/Av-k_ReImGammaL_ReEne_K_" << i + 1 <<  "_" << j + 1 << "_" << spin << ".dat";
                std::string var_left = oss_left.str();
                //std::cout << var_left << std::endl;
	            my_file_left.open(var_left, ios::in);
	    	    std::string line;

                int lineCount = -1;
                int energy_index = 0;
                while (std::getline(my_file_left, line)) {
                    lineCount++;
                    if (lineCount % parameters.grid_density == 0) {//this selects every 10th lines
                        if ((energy_index >= parameters.start.at(parameters.myid)) && (energy_index < parameters.start.at(parameters.myid) + parameters.steps_myid)) {
                            //this selects the appriotate lines for each process
                            std::istringstream iss(line);
                            double num1, num2, num3;
                            if (iss >> num1 >> num2 >> num3) {
                                // Use the values stored in num1, num2, and num3
                                int y = energy_index - parameters.start.at(parameters.myid);
                                gamma_left.at(y)(i, j) = num2 + parameters.j1 * num3;
                            } else {
                                std::cout << "Invalid line format: " << line << std::endl;
                            }
                        }

                        energy_index++;
                    }
                }
                my_file_left.close();

                fstream my_file_right;
                std::ostringstream oss_right;
                oss_right << parameters.path << "/Av-k_ReImGammaR_ReEne_K_" << i + 1 <<  "_" << j + 1 << "_" << spin << ".dat";
                std::string var_right = oss_right.str();
                //std::cout << var_left << std::endl;
	            my_file_right.open(var_right, ios::in);

                lineCount = -1;
                energy_index = 0;
                while (std::getline(my_file_right, line)) {
                    lineCount++;
                    if (lineCount % parameters.grid_density == 0) {//this selects every 10th lines
                        if ((energy_index >= parameters.start.at(parameters.myid)) && (energy_index < parameters.start.at(parameters.myid) + parameters.steps_myid)) {
                            //this selects the appriotate lines for each process
                            std::istringstream iss(line);
                            double num1, num2, num3;
                            if (iss >> num1 >> num2 >> num3) {
                                // Use the values stored in num1, num2, and num3
                                int y = energy_index - parameters.start.at(parameters.myid);
                                gamma_right.at(y)(i, j) = num2 + parameters.j1 * num3;
                            } else {
                                std::cout << "Invalid line format: " << line << std::endl;
                            }
                        }

                        energy_index++;
                    }
                }
                my_file_right.close();
            }
        }
    } else if (parameters.model_calc == 0) {
        for (int r = 0; r < parameters.steps_myid; r++) {
            for (int i = 0; i < parameters.num_orb_total; i++) {
                gamma_left.at(r)(i, i) = - 2.0 * parameters.gamma;
                gamma_right.at(r)(i, i) = - 2.0 * parameters.gamma;
            }
        }
    }   
}

void read_hamiltonian(const Parameters &parameters, MatrixType &hamiltonian_up, MatrixType &hamiltonian_down) {
    fstream my_file;
    std::ostringstream oss;
    oss << parameters.path << "/H.dat";
    std::string var = oss.str();
    //std::cout << var << std::endl;
	my_file.open(var, ios::in);
	std::string line;
    //cout << gf_non_int_r.at(i).size() << endl;
    if (my_file.is_open()) {
        std::string line;
        int lineCount = -1;
        while (std::getline(my_file, line)) {
            lineCount++;
            std::istringstream iss(line);
            double num1, num2, num3, num4, num5;
            if (iss >> num1 >> num2 >> num3 >> num4 >> num5) {
                // Use the values stored in num1, num2, and num3
                if (lineCount == 0) hamiltonian_up(0, 0) = num1;
                if (lineCount == 1) hamiltonian_down(0, 0) = num1;
            } 
        }
    } else {
        std::cout << var << " File cannot be opened for reading." << std::endl;
    }
    if (parameters.myid == 0) std::cout << hamiltonian_up << " " << hamiltonian_down << std::endl;
    my_file.close(); 
}

void get_interacting_retarded_gf(const Parameters &parameters, MatrixVectorType &gf_int_r, MatrixType &hamiltonian,
    const MatrixVectorType &sigma_mb_r, const MatrixVectorType &embedding_r, const MatrixVectorType &embedding_l) {
    if (parameters.read_gf == 1) {
        for (int r = 0; r < parameters.steps_myid; r++) {
            MatrixType inverse_gf = (gf_int_r.at(r).inverse() - sigma_mb_r.at(r));
            gf_int_r.at(r) = inverse_gf.inverse();
        }
        if (parameters.myid == 0) std::cout << "Got interacting green function from Dyson equation\n";
    } else {
        if (parameters.num_atoms != 1) {
            if (parameters.myid == 0) std::cout << "the number of atoms is more than 1. Can't construct gf with embedding self energies \n";
            exit(1);
        }
        Eigen::MatrixXcd gf_inverse(parameters.num_orbitals, parameters.num_orbitals);
        for (int r = 0; r < parameters.steps_myid; r++) {
            int y = parameters.start.at(parameters.myid) + r;
            gf_inverse(0, 0) = parameters.energy.at(y) - hamiltonian(0, 0) - embedding_r.at(r)(0, 0) - embedding_l.at(r)(0, 0) - sigma_mb_r.at(r)(0, 0);
            gf_int_r.at(r) = gf_inverse.inverse();
        }   
        if (parameters.myid == 0) std::cout << "Got interacting green function from the hamiltonian and self energies\n";
    }
}

void get_interacting_retarded_gf(const Parameters &parameters, MatrixVectorType &gf_int_r, MatrixType &hamiltonian,
    const MatrixVectorType &sigma_mb_r, const MatrixVectorType &delta) {
    Eigen::MatrixXcd gf_inverse(parameters.num_orbitals, parameters.num_orbitals);
    for (int r = 0; r < parameters.steps_myid; r++) {
        int y = parameters.start.at(parameters.myid) + r;
        gf_inverse(0, 0) = parameters.energy.at(y) - hamiltonian(0, 0) - delta.at(r)(0, 0);
        gf_int_r.at(r) = gf_inverse.inverse();
    }
}

void get_model_gf(const Parameters &parameters, MatrixVectorType &gf_int_r, const MatrixVectorType &sigma_mb_r,
    const MatrixVectorType &embedding_r, const MatrixVectorType &embedding_l) {
    for (int r = 0; r < parameters.steps_myid; r++) {
        int y = r + parameters.start.at(parameters.myid);
        gf_int_r.at(r)(0, 0) = 1.0 / (parameters.energy.at(y) + parameters.j1 * parameters.convergence - 
            parameters.onsite - sigma_mb_r.at(r)(0, 0) - embedding_l.at(r)(0, 0) - embedding_r.at(r)(0, 0));
    }
}

void get_lesser_se(Parameters &parameters, MatrixVectorType &se_lesser, MatrixVectorType &gamma, int left_right) {
    if (left_right == 0) {
        //cout << "getting a left self energy \n";
        for (int r = 0; r < parameters.steps_myid; r++) {
            int y = r + parameters.start.at(parameters.myid);
            se_lesser.at(r) = parameters.j1 * fermi_function(parameters.energy.at(y) - parameters.voltage, parameters) * gamma.at(r);
        }
    } else {
        //cout << "getting a right self energy \n";
        for (int r = 0; r < parameters.steps_myid; r++) {
            int y = r + parameters.start.at(parameters.myid);
            se_lesser.at(r) = parameters.j1 * fermi_function(parameters.energy.at(y) + parameters.voltage, parameters) * gamma.at(r);
        }
    }
}

void get_retarded_se(Parameters &parameters, MatrixVectorType &se_retarded, MatrixVectorType &gamma) {
    if (parameters.model_calc == 1) {//this is an ab initio calculation
        for (int r = 0; r < parameters.steps_myid; r++) {
            //double random = double((rand() % 1000)) / 1000 * 0;
            se_retarded.at(r)(0, 0) = - parameters.j1 * (gamma.at(r)(0, 0).real() / (2.0));
        }
    } else if (parameters.model_calc == 0) {
        std::vector<double> imag_self_energy_myid(parameters.steps_myid), imag_self_energy(parameters.steps);
        for (int r = 0; r < parameters.steps_myid; r++){
            imag_self_energy_myid.at(r) = - gamma.at(r)(0, 0).real() / (2.0);
        }

        MPI_Allgatherv(&(imag_self_energy_myid.at(0)), parameters.steps_myid, MPI_DOUBLE, &(imag_self_energy.at(0)),
	     &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE, MPI_COMM_WORLD);

	    for (int r = 0; r < parameters.steps_myid; r++) {
	    	int y = r + parameters.start.at(parameters.myid); 
	    	double impurity_self_energy_real_myid = kramer_kronig_relation(parameters, imag_self_energy, y);
	    	se_retarded.at(r)(0, 0) = impurity_self_energy_real_myid + parameters.j1 * imag_self_energy_myid.at(r);
	    }
    }
}

void get_lesser_gf(Parameters &parameters, MatrixVectorType &gf_int_r, MatrixVectorType &se_left,
     MatrixVectorType &se_right, MatrixVectorType &sigma_mb_l, MatrixVectorType &gf_int_l) {
    
    for(int r = 0; r < parameters.steps_myid; r++) {
        Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
        gf_int_l.at(r) = gf_int_r.at(r) * (se_left.at(r) + se_right.at(r) + sigma_mb_l.at(r)) * gf_int_a;
    } 
}