#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	//              omega   number of particles	    repulsion
	QuantumDots QD (1.0,	2, 						1);
	//				alpha0	alpha_step	alpha_max	beta0		beta_step	beta_max	system
	Investigate inv(0.5,   	0.1, 		4.0,  		0.0,	 	0.1, 		2.0,   		QD   	);
	//								 Monte carlo sim	Step length		Jastrow factor
	inv.compare_analytical_numerical(10000,			2.5,			1);
	inv.print_relative_energy_difference_to_file("energy_diff.csv");
	inv.print_alpha_meshgrid_to_file("test_alpha.csv");
	inv.print_beta_meshgrid_to_file("test_beta.csv");
}