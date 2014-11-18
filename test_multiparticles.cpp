#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(2);
	double alpha = 1.5, beta = 0.4, omega = 1;
	int number_of_particles = 2;
	int Jastrow = 1, Repulsion = 1;
	Trial_Wavefunction wf (alpha,beta,omega,number_of_particles,Jastrow);
	QuantumDots QD(omega,number_of_particles, Repulsion);
	QD.Set_Wavefunction(wf);

	vec expectation;

	//Number of Monte Carlo simulations
	int M = pow(10,6);
	
	cout << "M = " << M << endl;
	cout << "alpha = " << alpha << endl;
	cout << "beta = " << beta << endl;
	
	cout << "----------" << endl;
	cout << "Numerical: " << endl;
	expectation = QD.Importance_Sampling_Metropolis_Expectation_Values(M,0.2,1,0);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: " << endl;
	double tmp1 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl << endl; 

	cout << "----------" << endl;
	cout << "Analytical: " << endl;
	expectation = QD.Importance_Sampling_Metropolis_Expectation_Values(M,0.2,1,0);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: " << endl;
	double tmp2 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 
	


	cout << "------" << endl;
	cout << "Rel. Diff numerical-analytical: " << (tmp1-tmp2)/tmp1 << endl;
}