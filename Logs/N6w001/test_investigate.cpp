#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(12);
	// 				omega 	N 	R
	QuantumDots QD (0.01, 	6, 	1 );
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv(0.4, 	0.01, 		0.9,  	0.07, 	0.01, 	0.15,   		QD   );
	//				MCS	  			jastrow		ana. LE 	imp. samp. 		ana. QF 		delta_t
	inv.find_minimum(pow(10,5),		1, 			1,			0, 				0, 		 			0.1		);
	inv.print_energies_to_file("energies.csv");
	inv.print_variances_to_file("variances.csv");
	inv.print_alpha_meshgrid_to_file("alpha.csv");	
	inv.print_beta_meshgrid_to_file("beta.csv");
}
