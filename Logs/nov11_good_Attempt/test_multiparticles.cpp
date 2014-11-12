#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	double alpha = 1.15, beta = 0.09, omega = 1;
	int number_of_particles = 2;
	int Jastrow = 1, Repulsion = 1;
	Trial_Wavefunction wf (alpha,beta,omega,number_of_particles,Jastrow);
	QuantumDots QD(omega,number_of_particles, Repulsion);
	QD.Set_Wavefunction(wf);

	//Number of Monte Carlo simulations
	int M = 1000000;
	//Step length
	double delta_r = 2.5; //Gives around 50% acceptance rate for N=2
	
	cout << "M = " << M << endl;
	cout << "alpha = " << alpha << endl;
	cout << "beta = " << beta << endl;
	cout << "delta_r = " << delta_r << endl;
	
	cout << "Numerical: " << endl;
	vec expectation = QD.Brute_Force_Metropolis_Expectation_Values(M,delta_r,0);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: " << endl;
	double tmp1 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl << endl; 

	cout << "Analytical: " << endl;
	expectation = QD.Brute_Force_Metropolis_Expectation_Values(M,delta_r,1);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: " << endl;
	double tmp2 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 
	
	cout << "------" << endl;
	cout << "Rel. Diff numerical-analytical: " << (tmp1-tmp2)/tmp1 << endl;
}