#include "parameters.h"
#include <complex>  //this contains complex numbers and trig functions
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>
#include "utilis.h"
#include "gf.h"

void simple_tokenizer(std::string s, std::string &variable, std::string &value)
{
    std::stringstream ss(s);
    std::string word;
	int count = 0;
    while (ss >> word) {
        //std::cout << word << std::endl;
		if (count == 0) {
			variable = word;
		} else if (count == 1) {
			value = word;
		}
		count++;
    }
}


Parameters Parameters::from_file()
{
	Parameters parameters;
	std::string line, variable, value;
	std::ifstream input_file;

	input_file.open("input_file");

	if(!input_file.is_open())
	{
		std::cout << "The input_file doesn't exist \n" << std::endl;
		std::exit(1);
	} 

	parameters.se_mixing = 0.5;
	parameters.num_atoms = 1;
	parameters.read_se = 0;

	while (getline(input_file, line)) {
			//std::cout << line << '\n';
		simple_tokenizer(line, variable, value);
		//std::cout << "The variable name is " << variable << " The value is  " << value << std::endl;
		if (variable == "chemical_potential") {
			parameters.chemical_potential = std::stod(value);
		} else if (variable == "num_orbitals") {
			parameters.num_orbitals = std::stoi(value);
		} else if (variable == "voltage") {
			parameters.voltage = std::stod(value);
		} else if (variable == "interaction") {
			parameters.interaction = std::stoi(value);
		} else if (variable == "system") {
			parameters.system = std::stoi(value);
		} else if (variable == "path") {
			parameters.path = value;
		} else if (variable == "impurity_solver") {
			parameters.impurity_solver = std::stoi(value);
		} else if (variable == "hubbard_interaction") {
			parameters.hubbard_interaction = std::stod(value);
		} else if (variable == "self_consistent_steps") {
			parameters.self_consistent_steps = std::stoi(value);
		} else if (variable == "convergence") {
			parameters.convergence = std::stod(value);
		} else if (variable == "grid_density") {
			parameters.grid_density = std::stoi(value);
		} else if (variable == "read_gf") {
			parameters.read_gf = std::stoi(value);
		} else if (variable == "onsite") {
			parameters.onsite = std::stod(value);
		} else if (variable == "gamma") {
			parameters.gamma = std::stod(value);
		} else if (variable == "model_calc") {
			parameters.model_calc = std::stoi(value);
		}  else if (variable == "num_atoms") {
			parameters.num_atoms = std::stoi(value);
		} else if (variable == "se_mixing") {
			parameters.se_mixing = std::stod(value);
		} else if (variable == "read_se") {
			parameters.read_se = std::stoi(value);
		}	
		//else if (variable == "random_err") {
			//parameters.random_err = std::stod(value);
		//}
	}
	input_file.close();
	
	//parameters.energy.resize(parameters.steps);

	parameters.j1 = -1;
	parameters.j1 = sqrt(parameters.j1);

	//parameters.delta_energy = parameters.energy[1] - parameters.energy[0];
	//double delta_energy =
	//    (parameters.e_upper_bound - parameters.e_lower_bound) / (double)parameters.steps;
//
	//for (int i = 0; i < parameters.steps; i++) {
	//	parameters.energy[i] = parameters.e_lower_bound + delta_energy * (double)i;
	//}

	if (parameters.num_atoms != 1) parameters.read_gf = 1, parameters.model_calc = 1; //if there is more than 1 atom in the scattering region we must read the gf from file. 
	if (parameters.num_atoms != 1) parameters.read_se = 0; //if there is more than 1 atom in the scattering region we can't yet read the self energy. 
	

	if (parameters.read_gf == 1) parameters.model_calc = 1; //if we are reading the gf from file, then we can't do a model calc. 
	if (parameters.read_se == 1) parameters.impurity_solver = 1; //can't do kramer-kronig if we already have a self energy.


	parameters.voltage = 0.5 * parameters.voltage;
    MPI_Comm_size(MPI_COMM_WORLD, &parameters.size);
	MPI_Comm_rank(MPI_COMM_WORLD, &parameters.myid);

	parameters.energy.reserve(20100);

	get_number_energy_points(parameters);

	parameters.steps_proc.resize(parameters.size, 0), parameters.end.resize(parameters.size, 0), parameters.start.resize(parameters.size, 0),
		parameters.displs.resize(parameters.size, 0);

	decomp(parameters.steps, parameters.size, parameters.myid, &parameters.start[parameters.myid], &parameters.end[parameters.myid]);
	parameters.steps_myid = parameters.end[parameters.myid] - parameters.start[parameters.myid] + 1;

	parameters.steps_proc[parameters.myid] = parameters.steps_myid *  parameters.num_orbitals * parameters.num_orbitals;

	for (int a = 0; a < parameters.size; a++){
		MPI_Bcast(&parameters.start[a], 1, MPI_INT, a, MPI_COMM_WORLD);
		MPI_Bcast(&parameters.end[a], 1, MPI_INT, a, MPI_COMM_WORLD);
		MPI_Bcast(&parameters.steps_proc[a], 1, MPI_INT, a, MPI_COMM_WORLD);
	}

	parameters.num_orb_total = parameters.num_orbitals * parameters.num_atoms;

    parameters.displs[0] = 0;
    for (int i = 1; i < parameters.size; i++) {
        parameters.displs[i] = parameters.displs[i - 1] + parameters.steps_proc[i - 1];
    }
	

	if (parameters.myid == 0) {
		print_parameters(parameters);
	}

	parameters.delta_energy = parameters.energy[1] - parameters.energy[0];

	return parameters;
}

double fermi_function(double energy, const Parameters& parameters)
{
	return 1.0 / (1.0 + exp((energy - parameters.chemical_potential) / 0.001900069269));
}
//Parameters params = Parameters::from_file();

void print_parameters(Parameters& parameters)
{
	std::cout << "chemical_potential = " << parameters.chemical_potential << std::endl;
	//std::cout << "parameters.steps = " << parameters.steps << std::endl;
	std::cout << "num_orbitals = " << parameters.num_orbitals << std::endl;
	std::cout << "interaction is " << parameters.interaction << std::endl;
	std::cout << "system = " << parameters.system << std::endl;
	std::cout << "voltage = " << parameters.voltage << std::endl;
	std::cout << "steps = " << parameters.steps << std::endl;
	std::cout << "read_gf = " << parameters.read_gf << std::endl;
	std::cout << "model_calc = " << parameters.model_calc << std::endl;
	std::cout << "onsite = " << parameters.onsite << std::endl;
	std::cout << "gamma = " << parameters.gamma << std::endl;
	std::cout << "hubbard interaction = " << parameters.hubbard_interaction << std::endl;
	std::cout << "convergence = " << parameters.convergence << std::endl;
	std::cout << "max number of self consistent steps = " << parameters.self_consistent_steps << std::endl;
	std::cout << "grid density = " << parameters.grid_density << std::endl;
	std::cout << "num atoms = " << parameters.num_atoms << std::endl;
	std::cout << "num orb totals = " << parameters.num_orb_total << std::endl;
	std::cout << "se_mixing= " << parameters.se_mixing << std::endl;
	std::cout << "read_se= " << parameters.read_se << std::endl;

}
