#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	double alpha = 1, beta = 0, omega = 1;
	int number_of_particles = 6, 
	Trial_Wavefunction wf (alpha,beta,omega,number_of_particles,1);
	QuantumDots system(number_of_particles);
	system.Set_Hamiltonian(Unperturbed_Harmonic_Oscillator_Hamiltonian);
	system.Set_Wavefunction(wf);

	//Number of Monte Carlo simulations
	int M = 10;
	//Step length
	double delta_r = 2.5; //Gives around 50% acceptance rate
	cout << "delta_r = " << delta_r << endl;
	vec expect = Metropolis_Expectation_Values(system,M,delta_r);
	cout << expect;
	cout << "variance: " << abs(expect(1)-expect(0)*expect(0)) << endl;
}