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

	//Initialization
	double cumulative_function = 0; 		//Corresponds to local energy in the quantum example
	double cumulative_function_squared = 0;
	int counter = 0;

	while(counter < M){
		int i = rand() % 2; //Insert dimension here
		vec delta_vec_r = delta_r * randu<vec>(2);
		mat r_p = r;
		r_p.col(i) = r_p.col(i) + delta_vec_r;
		cout << "r:" << endl << r << endl;
		cout << "r_p:" << endl << r_p << endl;
		counter +=1;
	}

	//Create matrix for storing expectation values
	vec expectation_values = zeros(2);
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