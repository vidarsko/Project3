#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	//              omega   number of particles	    repulsion
	QuantumDots QD (1.0,	2, 						1);
	//				alpha0	alpha_step	alpha_max	beta0		beta_step	beta_max	system
	Investigate inv(0.86,   	0.01, 		0.96,  		-0.1,	 	0.01, 		0.1,   		QD   	);
	//				 Monte carlo sim	Step length		Jastrow factor
	inv.find_minimum(1000000,			2.5,			1);
	inv.print_energies_to_file("jas_rep_energies.csv");
	inv.print_alpha_meshgrid_to_file("jas_rep_alpha.csv");
	inv.print_beta_meshgrid_to_file("jas_rep_beta.csv");
}
