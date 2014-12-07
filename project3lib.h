#ifndef project2lib_h
#define project2lib_h
#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
#include <ctime>
using namespace std;
using namespace arma;



//*******************Trial Wavefunction class *****************//

class Trial_Wavefunction{
	/*
	The trial wavefunction object.
	*/
	public:
		//Trial parameters
		double alpha,beta;
		double omega;

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
		mat nabla_phi(double l, int k, mat r_i);
		double nabla2_phi(double l, int k, mat r_i);
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
		void Set_Wavefunction(Trial_Wavefunction wf);
		void update_omega(double w);

		//Local energy function
		double local_energy_function(mat r,int analytical_local_energy);

		//Help functions for the local energy
		double Potential(mat r);
		double Harmonic_potential(mat r);
		double Repulsive_potential(mat r);
		double Analytical_kinetic_energy(mat r);
		double numerical_sum_laplacians(Trial_Wavefunction wf,mat r,double h=1e-4);
		double LSP(Trial_Wavefunction wf, mat r, int i);
		double Average_distance(mat r);

		//Metropolis functions
		vec Brute_Force_Metropolis_Expectation_Values(int M, int analytical_local_energy);
		vec Importance_Sampling_Metropolis_Expectation_Values(int M, double delta_t, int analytical_local_energy, int analytical_quantum_force);
		vec Metropolis_interesting_quantities(int M);

		//Importance sampling help function
		mat quantum_force(mat r, int analytical_quantum_force, double h=1e-6);

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
		int alpha_dim, beta_dim;

		QuantumDots system;

		mat energies,variances, relative_energy_difference;

		//variables for the find parameters function
		vec alpha_optimal, beta_optimal, energy_optimal, variance_optimal, omegas;
		int number_of_omegas;

		//Variables for the analyze_importance_function
		mat imp_dt, imp_MCS, imp_energies, brute_energies; 

		//Variables for the compare_times function
		mat times, MCS_times;

		//Variables for the IQ-function
		double IQ_energy, IQ_variance, IQ_average_distance, IQ_potential_energy, IQ_harmonic_potential, IQ_repulsive_potential, IQ_kinetic_energy;


	public:

		//Constructor
		Investigate(double a0, double as, double am, double b0, double bs, double bm, QuantumDots sys);
		Investigate(QuantumDots sys);

		//Solve functions
		void find_minimum(int MCS, int jastrow, int analytical_local_energy,int importance_sampling, int analytical_quantum_force, double delta_t=0.1);
		void find_parameters(int MCS, vec omegas_input, int jastrow);
		void compare_analytical_numerical(int MCS, int jastrow);
		void analyze_importance_sampling(double alpha, double beta, int jastrow);
		void compare_times(double alpha, double beta, int jastrow);
		void interesting_quantities(double alpha, double beta);

		//Print functions
		void print_energies_to_file(string filename);
		void print_variances_to_file(string filename);
		void print_relative_energy_difference_to_file(string filename);
		void print_alpha_meshgrid_to_file(string filename);
		void print_beta_meshgrid_to_file(string filename);
		void print_optimals_to_file(string filename);
		void print_importance_analysis_to_files(string filenamebase);
		void print_times_to_file(string filenamebase);
		void print_interesting_quantities_to_file(string filename);
};









//**********Functions needed for class and elsewhere**********//


//Test functions
mat twobytwo(void);
mat twobysix(void);

#endif