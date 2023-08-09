
//#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cstdlib>
#include "parameters.h"


void get_spin_occupation(const Parameters &parameters, const std::vector<double> &gf_lesser_up,
                        const std::vector<double> &gf_lesser_down, double *spin_up, double *spin_down);

dcomp integrate(const Parameters& parameters, const std::vector<dcomp>& gf_1, const std::vector<dcomp>& gf_2, const std::vector<dcomp>& gf_3, const int r);

void self_energy_2nd_order(const Parameters& parameters, AIM &aim_up, AIM &aim_down);

double get_prefactor(const int i, const int j, const int r, const Parameters &parameters,
	 std::vector<double> &fermi_up, std::vector<double> &fermi_down);

double integrate_equilibrium(const Parameters& parameters, const std::vector<double>& gf_1, const std::vector<double>& gf_2,
    const std::vector<double>& gf_3, const int r, std::vector<double> &fermi_up, std::vector<double> &fermi_down);

void self_energy_2nd_order_kramers_kronig(const Parameters& parameters, AIM &aim_up, AIM &aim_down);

void impurity_solver_sigma_2(const Parameters &parameters, AIM &aim_up, AIM &aim_down);
