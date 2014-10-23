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
double two_particle_ground_state(mat r);
double laplacian_sum_ground_state(mat r);

//Laplacian functions
double sum_laplacian_i(double (*function),mat r,double);

//The local energy function 
double local_energy(double (*H)(double (*wf)(mat),mat),double (*wf)(mat),mat r);

//Hamiltonian functions
double Unperturbed_Harmonic_Oscillator_Hamiltonian(double (*wf)(mat), mat position);

//Hermite polynomials
double Hermite_polynomial(double x, int degree);



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



//*******************Trial Wavefunction class *****************//

class Trial_Wavefunction{
	/*
	The trial wavefunction object.
	*/
	private:
		double alpha,beta;
	public:

		//Constructor
		Trial_Wavefunction(double a, double b);
};
#endif