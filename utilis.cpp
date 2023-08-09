#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip> 
#include "parameters.h"
#include "utilis.h"
#include <mpi.h>
using namespace std;

void decomp(int steps, int size, int myid, int *s, int *e) {
    int remainder = steps % size; 
    int steps_per_proc = steps / size; //rounds towards 0
	if (myid < remainder) {
		*s = myid * (steps_per_proc + 1);
		*e = *s + steps_per_proc;
	} else {
		*s = myid * (steps_per_proc) + remainder;
		*e = *s + steps_per_proc - 1;
	}
}

void print_to_file(const Parameters &parameters, string filename, MatrixVectorType &quantity, int orbital_1, int orbital_2, int voltage_step){

	std::vector<dcomp> vec_1, vec_2;
	if (parameters.myid == 0) {
		//std::cout << "rank " << parameters.myid << " enters where parameters.myid == 0 \n";
		vec_1.resize(parameters.steps);
		vec_1.resize(parameters.steps);
		for (int r = 0; r < parameters.steps_myid; r ++){
			vec_1.at(r) = quantity.at(r)(orbital_1, orbital_2);  
		}
		for (int a = 1; a < parameters.size; a++){
			MPI_Recv(&vec_1.at(parameters.start.at(a)), parameters.steps_proc.at(a), MPI_DOUBLE_COMPLEX, a, 300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		//std::cout << "rank " << parameters.myid << " enters where parameters.myid != 0 \n";
		vec_2.resize(parameters.steps_myid);
		for (int r = 0; r < parameters.steps_myid; r++) {
			vec_2.at(r) = quantity.at(r)(orbital_1, orbital_2);
			//std::cout << "On rank " << parameters.myid << " vector has a value of " << vec_2.at(r)  << std::endl;
		}
		MPI_Send(&(vec_2.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, 0, 300, MPI_COMM_WORLD);
		//std::cout << "I, rank " << parameters.myid << " sent my part of the GF to rank 0 \n"; 
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (parameters.myid == 0) {

		std::ostringstream ossgf;
		ossgf  << "textfiles_from_code/" << voltage_step << "_" << filename << "_" << orbital_1 << "_" << orbital_2 << ".dat";
		std::string var = ossgf.str();
		std::ofstream file;
		//std::cout << "Printing file " << var <<  "\n";
	
		file.open(var);
		for (int r = 0; r < parameters.steps; r++) {
			file << parameters.energy.at(r) << "  " << vec_1.at(r).real() << "   " << vec_1.at(r).imag() << "\n";
		}
		file.close();					
	}
}

double kramer_kronig_relation(const Parameters &parameters, std::vector<double> &impurity_self_energy_imag, int r) {
    double se_real = 0;
    for (int i = 0; i < parameters.steps; i++) {
        if (i != r) {
            se_real += impurity_self_energy_imag.at(i) / (parameters.energy.at(i) - parameters.energy.at(r));
        }
    }
	return se_real * (parameters.energy.at(1) - parameters.energy.at(0)) / M_PI;
}

double absolute_value(double num1) {
	return std::sqrt((num1 ) * (num1));
}

void distribute_to_procs(const Parameters &parameters, std::vector<dcomp> &vec_1, std::vector<dcomp> &vec_2){
		MPI_Barrier(MPI_COMM_WORLD);

		if (parameters.myid == 0) {
			//std::cout << "rank " << parameters.myid << " enters where parameters.myid == 0 \n";
			for (int r = 0; r < parameters.steps_myid; r ++){
				vec_1.at(r) = vec_2.at(r);  
			}
			for (int a = 1; a < parameters.size; a++){
				MPI_Recv(&vec_1.at(parameters.start.at(a)), parameters.steps_proc.at(a), MPI_DOUBLE_COMPLEX, a, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//std::cout << "I, rank 0, recieved part of vec_1 from rank " << a << std::endl; 
				//for (int r = parameters.start.at(a); r < parameters.start.at(a) + parameters.steps_proc.at(a); r ++){
				//	std::cout << "This part has a value of " << vec_1.at(r) << " " << r << std::endl; 
				//}
			}
		} else {
			//std::cout << "rank " << parameters.myid << " enters where parameters.myid != 0 \n";
			MPI_Send(&(vec_2.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, 0, 200, MPI_COMM_WORLD);
			//std::cout << "I, rank " << parameters.myid << " sent my part of the GF to rank 0 \n"; 
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&(vec_1.at(0)), parameters.steps, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

}

void distribute_to_procs(const Parameters &parameters, std::vector<double> &vec_1, const std::vector<double> &vec_2){
		MPI_Barrier(MPI_COMM_WORLD);
		if (parameters.myid == 0) {
			//std::cout << "rank " << parameters.myid << " enters where parameters.myid == 0 \n";
			for (int r = 0; r < parameters.steps_myid; r ++){
				vec_1.at(r) = vec_2.at(r);  
			}
			for (int a = 1; a < parameters.size; a++){
				MPI_Recv(&vec_1.at(parameters.start.at(a)), parameters.steps_proc.at(a), MPI_DOUBLE, a, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//for (int r = parameters.start.at(a); r < parameters.start.at(a) + parameters.steps_proc.at(a); r ++){
				//	std::cout << "This part has a value of " << vec_1.at(r) << " " << r << std::endl; 
				//}
			}
		} else {
			//std::cout << "rank " << parameters.myid << " enters where parameters.myid != 0 \n";
			MPI_Send(&(vec_2.at(0)), parameters.steps_myid, MPI_DOUBLE, 0, 200, MPI_COMM_WORLD);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&(vec_1.at(0)), parameters.steps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

MatrixVectorType initializeMatrixVector(int size, int rows, int cols) {
    return MatrixVectorType(size, MatrixType::Zero(rows, cols));
}