#include <vector>
#include<iostream>
#include<fstream>
#include <string>
#include <complex>
#include "parameters.h"
using namespace std;
typedef std::complex<double> dcomp;

void split_c(string str, double* real_self_energy, double* imag_self_energy)
{
    string w = "";
    for (auto x : str)
    {
        if (x == ',')
        { 
            *real_self_energy = std::stod(w);
            w = "";
        }
        else {
            w = w + x;
        }
    }
    *imag_self_energy = std::stod(w);
}


string sanitize(string word)
{
int i = 0;

while(i < word.size())
{
    if(word[i] == '(' || word[i] == ')')
    {
        word.erase(i,1);
    } else {
        i++;
    }
}
return word;
}

void read_energy_grid(Parameters &parameters, string filename) {

    double real, imag;
    int number_of_lines = parameters.num_orbitals + 1;
    int number_orbitals_left = parameters.num_orbitals / 2.0;
	
    fstream my_file;

	std::ostringstream oss;
    oss << filename << "_reordered.dat";
    string var = oss.str();
	my_file.open(var, ios::in);

	if (!my_file) {
		cout << "No such file";
	} else {     
        int count = 0, energy_index = 0, orbital_index;
		string line;
		while (!my_file.eof()) {
			my_file >> line;
            orbital_index = count % (number_of_lines) - 1;
            energy_index = std::floor(count / number_of_lines);
            
            if (energy_index >= 2 * parameters.steps){
                break;
            }
            
            if (energy_index >= parameters.steps) {
                energy_index -= parameters.steps;
            }
            //std::cout << energy_index << std::endl;
            if (count % number_of_lines == 0){
                parameters.energy.at(energy_index) =  std::stod(line);
                //std::cout << "the energy.at(" << energy_index << ") is " << parameters.energy.at(energy_index) << std::endl;
            }  
            count++;
		}
	}
	my_file.close(); 


    //for (int r = 0; r < parameters.steps; r++) {
    //    std::cout << parameters.energy.at(r) << std::endl;
    //}

}


void read_self_energy(Parameters &parameters, std::vector<std::vector<dcomp>> &sigma_mb_up, std::vector<std::vector<dcomp>> &sigma_mb_down, string filename) {

    double real, imag;
    int number_of_lines = parameters.num_orbitals + 1;
    int number_orbitals_left = parameters.num_orbitals / 2.0;
	
    fstream my_file;

	std::ostringstream oss;
    oss << filename << "_reordered.dat";
    string var = oss.str();
	my_file.open(var, ios::in);

	if (!my_file) {
		cout << "No such file";
	} else {     
        int count = 0, energy_index = 0, orbital_index;
		string line;
		while (!my_file.eof()) {
			my_file >> line;
            orbital_index = count % (number_of_lines) - 1;
            energy_index = std::floor(count / number_of_lines);
            
            if (energy_index >= 2 * parameters.steps){
                break;
            }
            
            if (energy_index >= parameters.steps) {
                energy_index -= parameters.steps;
            }
            //std::cout << energy_index << std::endl;
            if (count % number_of_lines == 0){
                parameters.energy.at(energy_index) =  std::stod(line);
                //std::cout << "the energy.at(" << energy_index << ") is " << parameters.energy.at(energy_index) << std::endl;
            } else {
                //std::cout << orbital_index << " " << energy_index << "  " << sigma_mb_up.size() << "\n";
                line = sanitize(line);
                //std::cout << line << std::endl;
                split_c(line, &real, &imag);
                //std::cout << "The real part is still " << imag << std::endl;
                //std::cout << "Dividing by 30 gives " << imag/ (double)number_orbitals << std::endl;
                if (count > number_of_lines * parameters.steps) {
                    sigma_mb_down.at(orbital_index).at(energy_index) = real + parameters.j1 * imag;
                } else {
                    //cout << orbital_index << " " << energy_index << "\n";
                    sigma_mb_up.at(orbital_index).at(energy_index) = real + parameters.j1 * imag;
                    
                }         
            } 
            count++;
		}
	}
	my_file.close(); 
}



double reorder_self_energy(string file_name, const int number_orbitals){
    fstream my_file;
    std::ofstream myfile_reordered;

	std::ostringstream oss, oss1;
	oss << file_name << ".dat";
    oss1 << file_name << "_reordered.dat";
	std::string var = oss.str(), var1 = oss1.str();

    std::cout << "The files are " << var << " and " << var1 << std::endl;

	my_file.open(var, ios::in);
	myfile_reordered.open(var1);
    int count = 0, index = 0;
	if (!my_file) {
		cout << "No such file";
	} else {

		string line;
		while (!my_file.eof()) {
			my_file >> line;
            myfile_reordered << line << "\n";
            //std::cout << line << std::endl;
            count++;
            //std::cout << count << std::endl;
		}
	}

	my_file.close(); 
    myfile_reordered.close();
    //std::cout << (double)(count - 1) / ((number_orbitals + 1) * 2.0) << std::endl;
    return (double)(count - 1) / ((number_orbitals + 1) * 2.0);
}


void get_sigma(Parameters &parameters, std::vector<std::vector<dcomp>> &sigma_mb_r_up, std::vector<std::vector<dcomp>> &sigma_mb_r_down,
    std::vector<std::vector<dcomp>> &sigma_mb_l_up, std::vector<std::vector<dcomp>> &sigma_mb_l_down) {
    //g++  -O3 -g -std=c++11 average_self_energy.cpp -o average_self_energy 

    ifstream file("textfile_Fe_MgO/SigmaMBlesser.dat");
    bool sigma_lesser;
    if (!file) {// checks the existence of file
        std::cout << "SigmaMBlesser.dat file doesn't exist \n";
        sigma_lesser = false;
    } else {
		cout << "SigmaMBlesser.dat file exists \n";
        sigma_lesser = true;
    }

    parameters.steps = reorder_self_energy("textfile_Fe_MgO/SigmaMB", parameters.num_orbitals);

    if (sigma_lesser == true) {
        parameters.steps = reorder_self_energy("textfile_Fe_MgO/SigmaMBlesser", parameters.num_orbitals);
    }

    std::cout << "The number of energy points is " << parameters.steps << std::endl;

    parameters.energy.resize(parameters.steps);
    sigma_mb_r_up.resize(parameters.num_orbitals, std::vector<dcomp>(parameters.steps, 0));
    sigma_mb_r_down.resize(parameters.num_orbitals, std::vector<dcomp>(parameters.steps, 0));
    sigma_mb_l_up.resize(parameters.num_orbitals, std::vector<dcomp>(parameters.steps, 0));
    sigma_mb_l_down.resize(parameters.num_orbitals, std::vector<dcomp>(parameters.steps, 0));

    if (parameters.interaction == 1) {
        read_self_energy(parameters, sigma_mb_r_up, sigma_mb_r_down, "textfile_Fe_MgO/SigmaMB");
        read_self_energy(parameters, sigma_mb_l_up, sigma_mb_l_down, "textfile_Fe_MgO/SigmaMBlesser"); 
    } else {
        read_energy_grid(parameters, "textfile_Fe_MgO/SigmaMB"); //this si so the energy grid is created correctly.
    }
}



    //convert_to_electron_volts(self_energy_left_up, self_energy_right_up, self_energy_left_down, self_energy_right_down, number_energy_points);    
 