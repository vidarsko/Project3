#include "project3lib.h"


//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(int M, double delta_r, mat r){
	/*
	Function that uses the Monte Carlo Metropolis algorithm to 
	find the expectation value of some function g(r) with with 
	respect to a probabiltiy density P(r).
	Input:
		- int M 					- Number of Monte Carlo simulations
		- double delta_r 			- Predefined step length
		- mat r 					- Initial position r
		- function and prob.den 	- Not yet configured. 
	Output: 
		- vec expectation_values (2)- Expectation values of g and g^2. 
	*/

	cout << "Monte Carlo function properly initialized" << endl;
	cout << "Input parameters:" << endl;
	cout << "Number of monte carlo simulations M: " << M << endl;
	cout << "delta_r: " << delta_r << endl;
	cout << "Initial position r: " << endl << r << endl;
	vec expectation_values = zeros(2);
	return expectation_values;
}	

//Test functions for the metropolis algorithm 


//**********************Quantum Dots class *******************//

//Constructor
QuantumDots::QuantumDots(int a){
	N = a;
}

//Print functions 
void QuantumDots::print_N_to_terminal(void){
	cout << "N = " << N << endl;
}