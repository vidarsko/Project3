#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(20);

	vec omegas_input = zeros(6);
	omegas_input(0) = 0.01; omegas_input(1) = 0.10; omegas_input(2) = 0.28; 
	omegas_input(3) = 0.5; omegas_input(4) = 0.75; omegas_input(5) = 1;


	int MCS = pow(10,5);


	//N=2,  no jastrow
	cout << "N=2, no jastrow factor:" << endl;
	// 				omega 	N 	R
	QuantumDots QD1(1, 		2, 	1 );
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv1(0.01, 	0.01, 		1.2,  	0.01, 	0.01, 	0.03,   		QD1   );
	//					MCS			omegas_vector		jastrow
	inv1.find_parameters(MCS,		omegas_input, 		0);
	inv1.print_optimals_to_file("N2_nojastrow.csv");

	//N=2, jastrow
	cout << "N=2, with jastrow factor:" << endl;
	// 				omega 	N 	R
	QuantumDots QD2(1, 		2, 	1 );
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv2(0.01, 	0.01, 		1.2,  	0.01, 	0.01, 	1.0,   		QD2   );
	//					MCS			omegas_vector		jastrow
	inv2.find_parameters(MCS,		omegas_input, 		1);
	inv2.print_optimals_to_file("N2_jastrow.csv");		

	//N=6, no jastrow
	cout << "N=6, no jastrow factor:" << endl;
	// 				omega 	N 	R
	QuantumDots QD3(1, 		6, 	1 );
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv3(0.01, 	0.01, 		1.5,  	0.01, 	0.01, 	0.03,   		QD3   );
	//					MCS			omegas_vector		jastrow
	inv3.find_parameters(MCS,		omegas_input, 		0);
	inv3.print_optimals_to_file("N6_nojastrow.csv");	

	//N=6, jastrow
	cout << "N=6, with jastrow factor:" << endl;
	// 				omega 	N 	R
	QuantumDots QD4(1, 		6, 	1 );
	//				alpha0	alphastep	alpham 	beta0 	betas 	betamax 	system
	Investigate inv4(0.01, 	0.01, 		1.5,  	0.01, 	0.01, 	1,   		QD4   );
	//					MCS			omegas_vector		jastrow
	inv4.find_parameters(MCS,		omegas_input, 		1);
	inv4.print_optimals_to_file("N6_jastrow.csv");	
}