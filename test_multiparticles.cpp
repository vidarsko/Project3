#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	double alpha = 1, beta = 0, omega = 1;
	int number_of_particles = 6, spatial_dimension=2;
	Trial_Wavefunction wf (alpha,beta,omega,spatial_dimension,number_of_particles);
	QuantumDots system(number_of_particles);
	system.Set_Hamiltonian(Unperturbed_Harmonic_Oscillator_Hamiltonian);
	system.Set_Wavefunction(wf);

	mat r = randu<mat>(spatial_dimension,number_of_particles);
	//Number of Monte Carlo simulations
	int M = 10;
	//Step length
	double delta_r = 10; //Gives around 50% acceptance rate
	cout << "delta_r = " << delta_r << endl;
	vec expect = Metropolis_Expectation_Values(system,M,delta_r);
	cout << expect;
	cout << "variance: " << abs(expect(1)-expect(0)*expect(0)) << endl;
}