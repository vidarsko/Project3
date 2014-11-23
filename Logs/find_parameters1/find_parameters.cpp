#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(12);
	// 				omega 	N 	R
	QuantumDots QD (1, 	2, 	1 );
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv(0.0, 	0.1, 		1.2,  	0.0, 	0.01, 	1.2,   		QD   );

	vec omegas_input = zeros(3);
	omegas_input(0) = 0.01; omegas_input(1) = 0.28; omegas_input(2) = 1;

	int MCS = pow(10,5);
	//					MCS			omegas_vector	jastrow
	inv.find_parameters(MCS,	omegas_input, 	1);
	inv.print_optimals_to_file("jast_optimals.csv");
	inv.find_parameters(MCS, omegas_input, 0);
	inv.print_optimals_to_file("nojast_optimals.csv");


	//Time:
	/*
	alpha = 0: 	0.1:	1.2;
	beta = 0: 	0.1: 	1.2;
	MCS = pow(10,4)
	takes 2-3 minutes
	with pow(10,5): 7-8 minutes

	alpha = 0: 	0.1:	1.2;
	beta = 0: 	0.01: 	1.2;
	MCS = pow(10,5)
	takes 1hr 5 min

	*/
}
