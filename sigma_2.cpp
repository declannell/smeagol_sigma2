
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cstdlib>
#include "parameters.h"
#include "aim.h"
#include "utilis.h"
#include <iomanip>  



dcomp integrate(const Parameters& parameters, const std::vector<dcomp>& gf_1, const std::vector<dcomp>& gf_2, const std::vector<dcomp>& gf_3, const int r)
{

	dcomp result = 0;
	for (int i = 0; i < parameters.steps; i++) {
		for (int j = 0; j < parameters.steps; j++) {
			if (((i + j - r) > 0) && ((i + j - r) < parameters.steps)) {
				//this integrates the equation in PHYSICAL REVIEW B 74, 155125 2006
				//I say the green function is zero outside e_lower_bound and e_upper_bound. This means I need the final green function in the integral to be within an energy of e_lower_bound
				//and e_upper_bound. The index of 0 corresponds to e_lower_bound. Hence we need i+J-r>0 but in order to be less an energy of e_upper_bound we need i+j-r<steps.
				//These conditions ensure the enrgy of the gf3 greens function to be within (e_upper_bound, e_lower_bound)
				result += (parameters.delta_energy / (2.0 * M_PI)) * (parameters.delta_energy / (2.0 * M_PI)) * gf_1.at(i) * gf_2.at(j) * gf_3.at(i + j - r);
			}
		}
	}
	return result;
}

void self_energy_2nd_order(const Parameters& parameters, AIM &aim_up, AIM &aim_down)
{
	std::vector<dcomp> gf_retarded_up(parameters.steps), gf_retarded_down(parameters.steps), advanced_down(parameters.steps), 
		gf_lesser_up(parameters.steps), gf_lesser_down(parameters.steps), gf_greater_down(parameters.steps),
		advanced_down_myid(parameters.steps_myid), gf_lesser_up_myid(parameters.steps_myid), gf_lesser_down_myid(parameters.steps_myid),
		gf_greater_down_myid(parameters.steps_myid);

	for (int r = 0; r < parameters.steps_myid; r++) {
		advanced_down_myid.at(r) = std::conj(aim_down.impurity_gf_mb_retarded.at(r));
		gf_greater_down_myid.at(r) = aim_down.impurity_gf_mb_lesser.at(r) + aim_down.impurity_gf_mb_retarded.at(r) 
			- std::conj(aim_down.impurity_gf_mb_retarded.at(r));
	}

    MPI_Allgatherv(&(aim_up.impurity_gf_mb_retarded.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_retarded_up.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Allgatherv(&(aim_down.impurity_gf_mb_retarded.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_retarded_down.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Allgatherv(&(advanced_down_myid.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(advanced_down.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
	
    MPI_Allgatherv(&(aim_down.impurity_gf_mb_lesser.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_lesser_down.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Allgatherv(&(aim_up.impurity_gf_mb_lesser.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_lesser_up.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Allgatherv(&(gf_greater_down_myid.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_greater_down.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
	
	print_to_file(parameters,"gf_retarded_up", gf_retarded_up, 1);
	print_to_file(parameters,"advanced_down", advanced_down, 1);
	print_to_file(parameters,"gf_lesser_down", gf_lesser_down, 1);
	print_to_file(parameters,"gf_greater_down", gf_greater_down, 1);


	for (int r = 0; r < parameters.steps_myid; r++){
		int y = r + parameters.start.at(parameters.myid);
		//std::cout << "The energy index is " << y << std::endl;
		
		aim_up.self_energy_mb_retarded.at(r) = aim_up.self_energy_mb_retarded.at(r) * (1.0 - parameters.se_mixing) + parameters.se_mixing * parameters.hubbard_interaction * parameters.hubbard_interaction
		    * integrate(parameters, gf_retarded_up, gf_retarded_down, gf_lesser_down, y); //this resets the self energy	
		aim_up.self_energy_mb_retarded.at(r) += parameters.se_mixing * parameters.hubbard_interaction * parameters.hubbard_interaction
		    * integrate(parameters, gf_retarded_up, gf_lesser_down, gf_lesser_down, y); 
		aim_up.self_energy_mb_retarded.at(r) += parameters.se_mixing * parameters.hubbard_interaction * parameters.hubbard_interaction
		    * integrate(parameters, gf_lesser_up, gf_retarded_down, gf_lesser_down, y); 
		aim_up.self_energy_mb_retarded.at(r) += parameters.se_mixing * parameters.hubbard_interaction * parameters.hubbard_interaction
		    * integrate(parameters, gf_lesser_up, gf_lesser_down, advanced_down, y); 

		aim_up.self_energy_mb_lesser.at(r) = (1.0 - parameters.se_mixing) * aim_up.self_energy_mb_lesser.at(r) + parameters.se_mixing * (parameters.hubbard_interaction * parameters.hubbard_interaction
		    * integrate(parameters, gf_lesser_up, gf_lesser_down, gf_greater_down, y)); 	
	}
}

double get_prefactor(const int i, const int j, const int r, const Parameters &parameters,
	 std::vector<double> &fermi_up, std::vector<double> &fermi_down)
{
	int a = i + j - r;
	return fermi_up.at(i) * fermi_down.at(j) + fermi_down.at(a) * (1.0 - fermi_up.at(i) - fermi_down.at(j));
}

double integrate_equilibrium(const Parameters& parameters, const std::vector<double>& gf_1, const std::vector<double>& gf_2, 
	const std::vector<double>& gf_3, const int r, std::vector<double> &fermi_up, std::vector<double> &fermi_down)
{
	double result = 0;
	for (int i = 0; i < parameters.steps; i++) {
		for (int j = 0; j < parameters.steps; j++) {
			if (((i + j - r) > 0) && ((i + j - r) < parameters.steps)) {
				double prefactor = get_prefactor(i, j, r, parameters, fermi_up, fermi_down);
				//if (r == 123) count++;
				//this integrates the equation in PHYSICAL REVIEW B 74, 155125 2006
				//I say the green function is zero outside e_lower_bound and e_upper_bound. This means I need the final green function in the integral to be within an energy of e_lower_bound
				//and e_upper_bound. The index of 0 corresponds to e_lower_bound. Hence we need i+J-r>0 but in order to be less an energy of e_upper_bound we need i+j-r<steps.
				//These conditions ensure the enrgy of the gf3 greens function to be within (e_upper_bound, e_lower_bound)
				result += prefactor * (parameters.delta_energy / (M_PI)) * (parameters.delta_energy / (M_PI)) * gf_1.at(i) * gf_2.at(j) * gf_3.at(i + j - r);
			}
		}
	}
	return result;
}

//void self_energy_2nd_order_kramers_kronig(const Parameters& parameters, AIM &aim_up, AIM &aim_down)
//{
//	std::vector<double> impurity_self_energy_imag(parameters.steps); //this is for the kramer-kronig relation. 
//	std::vector<double> impurity_self_energy_real_myid(parameters.steps_myid), impurity_self_energy_imag_myid(parameters.steps_myid);
//
//	std::vector<double> fermi_eff_up(parameters.steps), fermi_eff_down(parameters.steps);
//	distribute_to_procs(parameters, fermi_eff_up,  aim_up.fermi_function_eff); //this is for the fermi function is the kk self energy.
//	distribute_to_procs(parameters, fermi_eff_down, aim_down.fermi_function_eff);
//
//	std::vector<dcomp> gf_lesser_up(parameters.steps), gf_lesser_down(parameters.steps), gf_greater_down(parameters.steps); //this is to pass into the integrater 
//		//for the self energy.
//	std::vector<dcomp> gf_lesser_up_myid(parameters.steps_myid), gf_lesser_down_myid(parameters.steps_myid), gf_greater_down_myid(parameters.steps_myid);
//	
//	for (int r = 0; r < parameters.steps_myid; r++) {
//		gf_lesser_down_myid.at(r) = aim_down.dynamical_field_lesser.at(r);
//		gf_lesser_up_myid.at(r) = aim_up.dynamical_field_lesser.at(r);
//		gf_greater_down_myid.at(r) = aim_down.dynamical_field_lesser.at(r) + aim_down.dynamical_field_retarded.at(r) - std::conj(aim_down.dynamical_field_retarded.at(r));
//	}	
//
//	distribute_to_procs(parameters, gf_lesser_up, gf_lesser_up_myid);
//	distribute_to_procs(parameters, gf_lesser_down, gf_lesser_down_myid);
//	distribute_to_procs(parameters, gf_greater_down, gf_greater_down_myid);	
//
//	std::vector<double> impurity_gf_up_imag(parameters.steps), impurity_gf_down_imag(parameters.steps);
//	std::vector<double> impurity_gf_up_imag_myid(parameters.steps_myid), impurity_gf_down_imag_myid(parameters.steps_myid);
//
//	for (int r = 0; r < parameters.steps_myid; r++) {
//		impurity_gf_up_imag_myid.at(r) = aim_up.dynamical_field_retarded.at(r).imag();
//		impurity_gf_down_imag_myid.at(r) = aim_down.dynamical_field_retarded.at(r).imag();
//    }
//
//	distribute_to_procs(parameters, impurity_gf_up_imag, impurity_gf_up_imag_myid);
//	distribute_to_procs(parameters, impurity_gf_down_imag, impurity_gf_down_imag_myid);	
//
//	for (int r = 0; r < parameters.steps_myid; r++){
//		int y = r + parameters.start.at(parameters.myid);
//		impurity_self_energy_imag_myid.at(r) = parameters.hubbard_interaction * parameters.hubbard_interaction
//		    * integrate_equilibrium(parameters, impurity_gf_up_imag, impurity_gf_down_imag, impurity_gf_down_imag, y, fermi_eff_up, fermi_eff_down); 	
//		
//		aim_up.self_energy_mb_lesser.at(r) = (parameters.hubbard_interaction * parameters.hubbard_interaction * (integrate(parameters, gf_lesser_up,
//			gf_lesser_down, gf_greater_down, y))); 
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//	distribute_to_procs(parameters, impurity_self_energy_imag, impurity_self_energy_imag_myid);
//
//	for (int r = 0; r < parameters.steps_myid; r++) {
//		int y = r + parameters.start.at(parameters.myid); 
//		impurity_self_energy_real_myid.at(r) = kramer_kronig_relation(parameters, impurity_self_energy_imag, y);
//		aim_up.self_energy_mb_retarded.at(r) = impurity_self_energy_real_myid.at(r) + parameters.j1 * impurity_self_energy_imag_myid.at(r);
//	}
//}
//
void self_energy_2nd_order_kramers_kronig(const Parameters& parameters, AIM &aim_up, AIM &aim_down)
{
	std::vector<double> impurity_self_energy_imag(parameters.steps); //this is for the kramer-kronig relation. 
	std::vector<double> impurity_self_energy_real_myid(parameters.steps_myid), impurity_self_energy_imag_myid(parameters.steps_myid);
	std::vector<double> fermi_eff_up(parameters.steps), fermi_eff_down(parameters.steps);
	std::vector<dcomp> gf_lesser_up(parameters.steps), gf_lesser_down(parameters.steps), gf_greater_down(parameters.steps); 
	std::vector<dcomp> gf_greater_down_myid(parameters.steps_myid);
	std::vector<double> impurity_gf_up_imag(parameters.steps), impurity_gf_down_imag(parameters.steps);
	std::vector<double> impurity_gf_up_imag_myid(parameters.steps_myid), impurity_gf_down_imag_myid(parameters.steps_myid);

	//MatrixVectorType test(parameters.steps_myid, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

	for (int r = 0; r < parameters.steps_myid; r++) {
		gf_greater_down_myid.at(r) =  aim_down.dynamical_field_lesser.at(r) 
			+ aim_down.dynamical_field_retarded.at(r) - std::conj(aim_down.dynamical_field_retarded.at(r));
		impurity_gf_up_imag_myid.at(r) = aim_up.dynamical_field_retarded.at(r).imag();
		impurity_gf_down_imag_myid.at(r) = aim_down.dynamical_field_retarded.at(r).imag();
	}	

    MPI_Allgatherv(&( aim_up.fermi_function_eff.at(0)), parameters.steps_myid, MPI_DOUBLE, &(fermi_eff_up.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgatherv(&( aim_down.fermi_function_eff.at(0)), parameters.steps_myid, MPI_DOUBLE, &(fermi_eff_down.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgatherv(&(aim_down.dynamical_field_lesser.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_lesser_down.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Allgatherv(&(aim_up.dynamical_field_lesser.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_lesser_up.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Allgatherv(&(gf_greater_down_myid.at(0)), parameters.steps_myid, MPI_DOUBLE_COMPLEX, &(gf_greater_down.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Allgatherv(&(impurity_gf_up_imag_myid.at(0)), parameters.steps_myid, MPI_DOUBLE, &(impurity_gf_up_imag.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgatherv(&(impurity_gf_down_imag_myid.at(0)), parameters.steps_myid, MPI_DOUBLE, &(impurity_gf_down_imag.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE, MPI_COMM_WORLD);

	//print_to_file(parameters, "gf_greater_down", gf_greater_down, 1);
	//print_to_file(parameters, "impurity_gf_up_imag", impurity_gf_up_imag, 1);
//
	//print_to_file(parameters, "gf_lesser_up", gf_lesser_up, 1);
//
	//if (parameters.myid == 0) {
	//	for (int r = 0; r < parameters.steps; r++) {
	//		std::cout << impurity_gf_up_imag.at(r) << std::endl;
	//	}
	//}


	for (int r = 0; r < parameters.steps_myid; r++){
		int y = r + parameters.start.at(parameters.myid);
		impurity_self_energy_imag_myid.at(r) = parameters.hubbard_interaction * parameters.hubbard_interaction
		    * integrate_equilibrium(parameters, impurity_gf_up_imag, impurity_gf_down_imag, impurity_gf_down_imag, y, fermi_eff_up, fermi_eff_down); 	
		
		aim_up.self_energy_mb_lesser.at(r) = (parameters.hubbard_interaction * parameters.hubbard_interaction * (integrate(parameters, gf_lesser_up,
			gf_lesser_down, gf_greater_down, y))); 
		//test.at(r)(0, 0) = impurity_self_energy_imag_myid.at(r);
	}	

	//print_to_file(parameters, "test", test, 0, 0, 1);
		//std::cout << aim_up.dynamical_field_retarded.at(r) << " " << impurity_gf_up_imag.at(r) << " " << parameters.energy.at(r) << std::endl;

    MPI_Allgatherv(&(impurity_self_energy_imag_myid.at(0)), parameters.steps_myid, MPI_DOUBLE, &(impurity_self_energy_imag.at(0)),
	 &(parameters.steps_proc.at(0)), &(parameters.displs.at(0)), MPI_DOUBLE, MPI_COMM_WORLD);

	for (int r = 0; r < parameters.steps_myid; r++) {
		int y = r + parameters.start.at(parameters.myid); 
		impurity_self_energy_real_myid.at(r) = kramer_kronig_relation(parameters, impurity_self_energy_imag, y);
		aim_up.self_energy_mb_retarded.at(r) = impurity_self_energy_real_myid.at(r) + parameters.j1 * impurity_self_energy_imag_myid.at(r);
	}
}


void get_difference_self_energy(const Parameters &parameters, std::vector<dcomp> &self_energy_mb_up,
	std::vector<dcomp> &old_self_energy_mb_up, double &difference){
	double difference_proc = - std::numeric_limits<double>::infinity();
	double real_difference = 0, imag_difference = 0;
	for (int r = 0; r < parameters.steps_myid; r++) {
		real_difference = absolute_value(self_energy_mb_up.at(r).real() - old_self_energy_mb_up.at(r).real());
		imag_difference = absolute_value(self_energy_mb_up.at(r).imag() - old_self_energy_mb_up.at(r).imag());
		//std::cout << gf_local_up.at(r)(i, j).real() << " " << old_green_function.at(r)(i, j).real() << std::endl;
		//std::cout << real_difference << "  " << imag_difference << "  "  << difference << "\n";
		difference_proc = std::max(difference_proc, std::max(real_difference, imag_difference));
		old_self_energy_mb_up.at(r) = self_energy_mb_up.at(r);
	}
	//std::cout << "I am rank " << parameters.myid << ". The difference for me is " << difference_proc << std::endl;
	//MPI_Allreduce would do the same thing.
	MPI_Allreduce(&difference_proc, &difference, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
}




void impurity_solver_sigma_2(const Parameters &parameters, AIM &aim_up, AIM &aim_down) {

	if (parameters.impurity_solver == 2) {// kramer kronig relation.
		if (parameters.myid == 0) std::cout << "using the kramer-kronig relation for second order perturbation theory\n";
		self_energy_2nd_order_kramers_kronig(parameters, aim_up, aim_down);
		self_energy_2nd_order_kramers_kronig(parameters, aim_down, aim_up);
	} else if (parameters.impurity_solver == 1) {//brute force sigma_2
		if (parameters.myid == 0) std::cout << "using the brute force relation for second order perturbation theory\n";
		double difference = std::numeric_limits<double>::infinity();
		int count = 0;
		std::vector<dcomp> old_self_energy_mb(parameters.steps_myid, 0);
		while (difference > parameters.convergence && count < parameters.self_consistent_steps) {
			self_energy_2nd_order(parameters, aim_up, aim_down);
			self_energy_2nd_order(parameters, aim_down, aim_up);
			
			get_difference_self_energy(parameters, aim_up.self_energy_mb_retarded, old_self_energy_mb, difference);
			if (parameters.myid == 0) std::cout << std::setprecision(15) << "The difference is " << difference << ". The count is " << count << std::endl;

			aim_up.get_impurity_gf_mb_retarded(parameters);
			aim_down.get_impurity_gf_mb_retarded(parameters);
			aim_up.get_impurity_gf_mb_lesser(parameters);
			aim_down.get_impurity_gf_mb_lesser(parameters);

			if (difference < parameters.convergence) break;
			count++;
		}
	}
}
