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


//*******************Trial Wavefunction class *****************//

class Trial_Wavefunction{
	/*
	The trial wavefunction object.
	*/
	private:
		//Trial parameters
		double alpha,beta;
		double omega;

	public:
		//System settings
		int number_of_particles;
		int spatial_dimension;

		//Constructor
		Trial_Wavefunction(double a, double b, double c, int S,int N);
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
	public:
		int number_of_particles;  //Number of electrons in trap
		double (*Hamiltonian) (Trial_Wavefunction wf, mat r); //Hamiltonian
		Trial_Wavefunction Wave_function;

	
		//Constructor
		QuantumDots(int N);

		//Uncategorized
		void Set_Hamiltonian(double (*H)(Trial_Wavefunction wf, mat r));
		void Set_Wavefunction(Trial_Wavefunction wf);
		double local_energy(mat r);

		//Print functions 
		void print_numberofparticles_to_terminal(void);
};

//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(QuantumDots system, int M, double delta_r);

//Test functions for the metropolis algorithm 
double Test_Probability_Density(mat r);
double Test_Evaluation_Function(mat r);
double two_particle_ground_state(mat r);

//Laplacian functions
double sum_laplacians(Trial_Wavefunction wf,mat r,double h=1e-4);

//Hamiltonian functions
double Unperturbed_Harmonic_Oscillator_Hamiltonian(Trial_Wavefunction Wave_function, mat r);


#endif