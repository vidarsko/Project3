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

//Hamiltonian functions
double Unperturbed_Harmonic_Oscillator_Hamiltonian(double (*wf)(mat), mat r);

//*******************Trial Wavefunction class *****************//

class Trial_Wavefunction{
	/*
	The trial wavefunction object.
	*/
	private:
		//Trial parameters
		double alpha,beta;

		//System settings
		double omega;
		int number_of_particles;
		int spatial_dimension;

	public:
		//Constructor
		Trial_Wavefunction(double a, double b,double c,int S,int N);
		Trial_Wavefunction();

		//Call functions 
		double call(mat r);
		double call_squared(mat r);

		//Help functions
		double Hermite_polynomial(double x, int degree);
		double nx(int i);
		double ny(int i);
		double phi(int i,mat r_i);

};



//**********************Quantum Dots class *******************//

class QuantumDots{
	/*
	The system
	*/
	private:
		int number_of_particles;  //Number of electrons in trap
		double Hamiltonian (double (*wf)(mat), mat r); //Hamiltonian
		Trial_Wavefunction Wave_function;

	public:
		//Constructor
		QuantumDots(int N);

		//Uncategorized
		void Set_Hamiltonian(double (*H) (double (*wf)(mat), mat r));
		void Set_Wavefunction(Trial_Wavefunction wf);
		double Wave_function_eval(mat r);
		double local_energy(mat r);
		double laplacian_sum_ground_state(mat r);

		//Print functions 
		void print_numberofparticles_to_terminal(void);
};

#endif