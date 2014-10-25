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
		int Jastrow_factor;

		//Constructor
		Trial_Wavefunction(double a, double b, double c, int N, int J);
		Trial_Wavefunction();

		//Call functions 
		double call(mat r);
		double call_squared(mat r);

		//Help functions
		double Hermite_polynomial(double x, int degree);
		double nx(int i);
		double ny(int i);
		double phi(int i,mat r_i);
		double a(int i, int j);

};



//**********************Quantum Dots class *******************//

class QuantumDots{
	/*
	The system
	*/
	public:
		double omega;
		int number_of_particles;  //Number of electrons in trap
		int repulsion; //Repulsion part of hamiltonian?
		Trial_Wavefunction Wave_function;

	
		//Constructor
		QuantumDots(double w, int N, int R);
		QuantumDots();

		//Configuration functions
		double Hamiltonian (Trial_Wavefunction wf, mat r); //Hamiltonian
		void Set_Wavefunction(Trial_Wavefunction wf);

		//Local energy function
		double local_energy(mat r);

		//Print functions 
		void print_numberofparticles_to_terminal(void);
};



//********************Investigation class **************************//

class Investigate{
	/*
	Class made to find the optimal parameters when determining the expectation value of the energy.
	*/
	private:
		double alpha_0,alpha_step,alpha_max;
		double beta_0,beta_step,beta_max;

		QuantumDots system;

		mat energies, variances;
	public:

		//Constructor
		Investigate(double aa, double as, double am, double b0, double bs, double bm, QuantumDots sys);

		//Solve function
		void solve(int MCS, double delta_r, int jastrow);

		//Print functions
		void print_energies_to_file(string filename);
		void print_variances_to_file(string filename);
		void print_alpha_meshgrid_to_file(string filename);
		void print_beta_meshgrid_to_file(string filename);

};









//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(QuantumDots system, int M, double delta_r);

//Laplacian functions
double sum_laplacians(Trial_Wavefunction wf,mat r,double h=1e-4);



#endif