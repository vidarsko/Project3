#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(20);

	vec omegas_input = zeros(3);
	omegas_input(0) = 0.01; omegas_input(1) = 0.28; omegas_input(2) = 1.0; 

	int MCS = pow(10,5);

	//N=6, jastrow
	cout << "N=12, with jastrow factor:" << endl;
	// 				omega 	N 	R
	QuantumDots QD4(1, 		12, 	1 );
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv4(0.01, 	0.01, 		1.5,  	0.01, 	0.01, 	1,   		QD4   );
	//					MCS			omegas_vector		jastrow
	inv4.find_parameters(MCS,		omegas_input, 		1);
	inv4.print_optimals_to_file("N12_jastrow.csv");	
}
