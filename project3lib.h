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
	private:
		double omega;

	public:
		int number_of_particles;  //Number of electrons in trap
		Trial_Wavefunction Wave_function;

	
		//Constructor
		QuantumDots(double a, int N);

		//Configuration functions
		double Hamiltonian (Trial_Wavefunction wf, mat r, int R); //Hamiltonian
		void Set_Wavefunction(Trial_Wavefunction wf);

		//Local energy function
		double local_energy(mat r, int R);

		//Print functions 
		void print_numberofparticles_to_terminal(void);
};

//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(QuantumDots system, int M, double delta_r, int R=0);

//Laplacian functions
double sum_laplacians(Trial_Wavefunction wf,mat r,double h=1e-4);



#endif