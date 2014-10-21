#include "project3lib.h"


//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(double (*P)(mat), double (*g)(mat), int M, double delta_r, mat r){
	/*
	Function that uses the Monte Carlo Metropolis algorithm to 
	find the expectation value of some function g(r) with with 
	respect to a probabiltiy density P(r).
	Input:
		- double P 					- Funciton propto Prob density P
		- double g					- Function whose exp. value to be evaluated.
		- int M 					- Number of Monte Carlo simulations
		- double delta_r 			- Predefined step length
		- mat r 					- Initial position r
	Output: 
		- vec expectation_values (2)- Expectation values of g and g^2. 
	*/

	/*
	cout << "Monte Carlo function properly initialized" << endl;
	cout << "Input parameters:" << endl;
	cout << "Number of monte carlo simulations M: " << M << endl;
	cout << "delta_r: " << delta_r << endl;
	cout << "Initial position r: " << endl << r << endl;
	*/
	vec expectation_values = zeros(2);
	expectation_values(0) = Test_Evaluation_Function(r);
	expectation_values(1) = pow(Test_Evaluation_Function(r),2);
	cout << "Expectation values" << endl << expectation_values << endl;
	return expectation_values;
}	

//Test functions for the metropolis algorithm 
double Test_Probability_Density(mat r){
	/*
	Probability density where P(r)=e^(-radius) where radius = |mat r|:
	*/
	double radius = norm(r,2);
	return exp(-radius);
	
}

double Test_Evaluation_Function(mat r){
	/*
	Function to be evaluated in the probability density. 
	Here, the value is 
	*/
	return norm(r,2);
}

















//**********************Quantum Dots class *******************//

//Constructor
QuantumDots::QuantumDots(int a){
	N = a;
}

//Print functions 
void QuantumDots::print_N_to_terminal(void){
	cout << "N = " << N << endl;
}