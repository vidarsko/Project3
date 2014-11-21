#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	// 				omega 	N 	R
	QuantumDots QD (1.5, 	6, 	0);
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv(0.9, 	0.05, 		1.2,  	0.0, 	0.1, 	0.1,   		QD   );
	//				MCS	  			jastrow		ana. LE 	imp. samp. 		ana. QF 		delta_t
	inv.find_minimum(pow(10,7),		0, 			0,			0, 				0, 		 			0.1		);
	inv.print_energies_to_file("energies.csv");
	inv.print_variances_to_file("variances.csv");
	inv.print_alpha_meshgrid_to_file("alpha.csv");
	inv.print_beta_meshgrid_to_file("beta.csv");
}
