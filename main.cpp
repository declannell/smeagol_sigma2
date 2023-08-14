#include "parameters.h"
#include "gf.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip>  
#include <string.h>
#include "read_sigma.h"
#include </usr/include/eigen3/Eigen/Dense>
#include "current.h"
#include "utilis.h"
#include <mpi.h>
#include "aim.h"
#include "sigma_2.h"


int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

	Parameters parameters = Parameters::from_file();
    if (parameters.myid == 0) std::cout << std::setprecision(15) << "created parameters \n";

	int num_orb_total = parameters.num_orb_total, steps_myid = parameters.steps_myid;

    MatrixVectorType gf_int_r_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType gf_int_r_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);

    MatrixVectorType sigma_mb_r_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType sigma_mb_r_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType sigma_mb_l_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType sigma_mb_l_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);

    MatrixVectorType gamma_left_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType gamma_right_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType gamma_left_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType gamma_right_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);

    MatrixVectorType se_left_lesser_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType se_right_lesser_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType se_left_lesser_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType se_right_lesser_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);

    MatrixVectorType se_left_retarded_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType se_right_retarded_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType se_left_retarded_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType se_right_retarded_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);

    MatrixVectorType gf_int_l_up = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);
    MatrixVectorType gf_int_l_down = initializeMatrixVector(steps_myid, num_orb_total, num_orb_total);

    MatrixType hamiltonian_up(num_orb_total, num_orb_total);
    MatrixType hamiltonian_down(num_orb_total, num_orb_total);

	if (parameters.myid == 0) std::cout << "initialised matrix \n" ;
	read_gamma(parameters, gamma_left_up, gamma_right_up, 1);
	read_gamma(parameters, gamma_left_down, gamma_right_down, 2);

	if (parameters.myid == 0) std::cout << "Got the gamma matrices from file \n";

	get_lesser_se(parameters, se_left_lesser_up, gamma_left_up, 0);
	get_lesser_se(parameters, se_left_lesser_down, gamma_left_down, 0);
	get_lesser_se(parameters, se_right_lesser_up, gamma_right_up, 1);
	get_lesser_se(parameters, se_right_lesser_down, gamma_right_down, 1);

	if (parameters.myid == 0) cout << "Got the lesser self energy using FD \n";

	if (parameters.read_gf == 1) {//if we are reading the green function from file we don't use the retarded self energy. 
		read_non_interacting_gf(parameters, gf_int_r_up, 1);
		read_non_interacting_gf(parameters, gf_int_r_down, 2);
		if (parameters.myid == 0) std::cout << "Read the non-interacting green function from file \n";
	} else {//can either do a model calculation or we are creating it from the hamiltonian and embedding self energies
		//if we don't read the gf from file we need the retarded embedding self energies
		get_retarded_se(parameters, se_left_retarded_up, gamma_left_up);
		get_retarded_se(parameters, se_left_retarded_down, gamma_left_down);
		get_retarded_se(parameters, se_right_retarded_up, gamma_right_up);
		get_retarded_se(parameters, se_right_retarded_down, gamma_right_down);

		if (parameters.myid == 0) cout << "Got the retarded self energy from the kramer kronig \n";

		if (parameters.model_calc == 0) {
			if (parameters.myid == 0) std::cout << "Doing a model calculation \n";
			get_model_gf(parameters, gf_int_r_up, sigma_mb_r_up, se_right_retarded_up, se_left_retarded_up);
    		get_model_gf(parameters, gf_int_r_down, sigma_mb_r_down, se_right_retarded_down, se_left_retarded_down);
		} else if (parameters.model_calc == 1) {
			read_hamiltonian(parameters, hamiltonian_up, hamiltonian_down);
			if (parameters.myid == 0) std::cout << "Read the hamiltonian form file \n";
    		get_interacting_retarded_gf(parameters, gf_int_r_up, hamiltonian_up, sigma_mb_r_up, se_right_retarded_up, se_left_retarded_up);
    		get_interacting_retarded_gf(parameters, gf_int_r_down, hamiltonian_down, sigma_mb_r_down, se_right_retarded_down, se_left_retarded_down);
		}
	}

	for (int i = 0; i < parameters.num_orb_total; i++) {
		for (int j = 0; j < parameters.num_orb_total; j++) {
			print_to_file(parameters, "gamma_left_up", gamma_left_up, i, j , 1);
			print_to_file(parameters, "gamma_right_up", gamma_right_up, i, j , 1);
			print_to_file(parameters, "gamma_left_down", gamma_left_down, i, j , 1);
			print_to_file(parameters, "gamma_right_down", gamma_right_down, i, j , 1);
			if (parameters.read_gf != 1) {
				print_to_file(parameters, "se_left_retarded_up", se_left_retarded_up, i, j , 1);
				print_to_file(parameters, "se_left_retarded_down", se_left_retarded_down, i, j , 1);
				print_to_file(parameters, "se_right_retarded_up", se_right_retarded_up, i, j , 1);
				print_to_file(parameters, "se_right_retarded_down", se_right_retarded_down, i, j , 1);
			}
			print_to_file(parameters, "gf_non_r_up", gf_int_r_up, i, j , 1);
			print_to_file(parameters, "gf_non_r_down", gf_int_r_down, i, j , 1);
		}
	}

	get_lesser_gf(parameters, gf_int_r_up, se_left_lesser_up, se_right_lesser_up, sigma_mb_l_up, gf_int_l_up);
	get_lesser_gf(parameters, gf_int_r_down, se_left_lesser_down, se_right_lesser_down, sigma_mb_l_down, gf_int_l_down);

	for (int i = 0; i < parameters.num_atoms; i++) {
		AIM aim_up(parameters, gf_int_r_up, gf_int_l_up, i);
		AIM aim_down(parameters, gf_int_r_down, gf_int_l_down, i);

		if (parameters.hubbard_interaction != 0) {
			impurity_solver_sigma_2(parameters, aim_up, aim_down);
		}

		for (int r = 0; r < parameters.steps_myid; r++) {
			sigma_mb_r_up.at(r) (i, i) = aim_up.self_energy_mb_retarded.at(r);
			sigma_mb_l_up.at(r) (i, i) = aim_up.self_energy_mb_lesser.at(r);
			sigma_mb_r_down.at(r) (i, i) = aim_down.self_energy_mb_retarded.at(r);
			sigma_mb_l_down.at(r) (i, i) = aim_down.self_energy_mb_lesser.at(r);
		}
	}

	if (parameters.model_calc == 0) {//doing a model calculation
		get_model_gf(parameters, gf_int_r_up, sigma_mb_r_up, se_right_retarded_up, se_left_retarded_up);
    	get_model_gf(parameters, gf_int_r_down, sigma_mb_r_down, se_right_retarded_down, se_left_retarded_down);
	} else if (parameters.model_calc == 1) {//doing an ab initio calculation
		get_interacting_retarded_gf(parameters, gf_int_r_up, hamiltonian_up, sigma_mb_r_up, se_right_retarded_up, se_left_retarded_up);
    	get_interacting_retarded_gf(parameters, gf_int_r_down, hamiltonian_down, sigma_mb_r_down, se_right_retarded_down, se_left_retarded_down);
	}

	get_lesser_gf(parameters, gf_int_r_up, se_left_lesser_up, se_right_lesser_up, sigma_mb_l_up, gf_int_l_up);
	get_lesser_gf(parameters, gf_int_r_down, se_left_lesser_down, se_right_lesser_down, sigma_mb_l_down, gf_int_l_down);

	for (int i = 0; i < parameters.num_orb_total; i++) {
		print_to_file(parameters, "se_retarded_up", sigma_mb_r_up, i, i , 1);
		print_to_file(parameters, "se_lesser_up", sigma_mb_l_up, i, i , 1);
		print_to_file(parameters, "se_retarded_down", sigma_mb_r_down, i, i , 1);
		print_to_file(parameters, "se_lesser_down", sigma_mb_l_down, i, i , 1);
		for (int j = 0; j < parameters.num_orb_total; j++) {
			print_to_file(parameters, "gf_retarded_up", gf_int_r_up, i, j , 1);
			print_to_file(parameters, "gf_lesser_up", gf_int_l_up, i, j , 1);
			print_to_file(parameters, "gf_retarded_down", gf_int_r_down, i, j , 1);
			print_to_file(parameters, "gf_lesser_down", gf_int_l_down, i, j , 1);
		}
	}

	double current_up_l, current_up_r, current_down_l, current_down_r;
	double current_up_coherent, current_down_coherent;

	get_current(parameters, gf_int_r_up, gf_int_l_up, gamma_left_up, current_up_l, 0);
	get_current(parameters, gf_int_r_up, gf_int_l_up, gamma_right_up, current_up_r, 1);
	get_current(parameters, gf_int_r_down, gf_int_l_down, gamma_left_down, current_down_l, 0);
	get_current(parameters, gf_int_r_down, gf_int_l_down, gamma_right_down, current_down_r, 1);

	get_current_transmission(parameters, gf_int_r_up, gamma_left_up, gamma_right_up, current_up_coherent);
	get_current_transmission(parameters, gf_int_r_down, gamma_left_down, gamma_right_down, current_down_coherent);

	if (parameters.myid == 0) {
		cout << std::setprecision(15) << "The spin up Meir wingreen left current is " << current_up_l << endl;
		cout << "The spin up Meir wingreen right current is " << current_up_r << endl;
		cout << "The coherent spin up current is " << current_up_coherent << endl;

		cout << "The spin down Meir wingreen left current is " << current_down_l << endl;
		cout << "The spin down Meir wingreen right current is " << current_down_r << endl;
		cout << "The coherent spin down current is " << current_down_coherent << endl;
	}

	MPI_Finalize();
}