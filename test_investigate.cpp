#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	// 				omega 	N 	R
	QuantumDots QD (1.0, 	2, 	1);
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv(0.7, 	0.1, 		1.2,  	0.2, 	0.01, 	0.5,   		QD   );
	//				MCS	  		jastrow		ana. LE 	imp. samp. 		ana. QF 		delta_t
	inv.find_minimum(10000,		1, 			0,			1, 				1, 		 			0.1		);
	inv.print_energies_to_file("test_energies.csv");
	inv.print_variances_to_file("test_variances.csv");
	inv.print_alpha_meshgrid_to_file("test_alpha.csv");
	inv.print_beta_meshgrid_to_file("test_beta.csv");
}
