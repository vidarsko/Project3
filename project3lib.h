#ifndef project2lib_h
#define project2lib_h
#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
using namespace std;
using namespace arma;

//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(double (*P)(mat), double (*g)(mat), int M, double delta_r, mat r);


//Test functions for the metropolis algorithm 
double Test_Probability_Density(mat r);
double Test_Evaluation_Function(mat r);

//**********************Quantum Dots class *******************//

class QuantumDots{
	/*
	The system.
	*/
	private:
		int N; //Number of electrons in trap.
	public:

		//Constructor
		QuantumDots(int a);

		//Print functions 
		void print_N_to_terminal(void);

};

#endif