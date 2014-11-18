#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(1);
	double alpha = 1.01, beta = 0.4, omega = 1;
	int number_of_particles = 2;
	int Jastrow = 0, Repulsion = 0;
	double delta_t = 1e-4;
	Trial_Wavefunction wf (alpha,beta,omega,number_of_particles,Jastrow);
	QuantumDots QD(omega,number_of_particles, Repulsion);
	QD.Set_Wavefunction(wf);

	vec expectation;

	//Number of Monte Carlo simulations
	int M = pow(10,6);
	
	cout << "M = " << M << endl;
	cout << "alpha = " << alpha << endl;
	cout << "beta = " << beta << endl;
	cout << "delta_t = " << delta_t << endl;
	double tmp1, tmp2;
	/*
	cout << "----------" << endl;
	cout << "Importance sampling, numerical local energy, numerical quantum force: " << endl;
	expectation = QD.Importance_Sampling_Metropolis_Expectation_Values(M,delta_t,0,0);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: ";
	tmp1 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 

	cout << "----------" << endl;
	cout << "Importance sampling, numerical local energy, analytical quantum force: " << endl;
	expectation = QD.Importance_Sampling_Metropolis_Expectation_Values(M,delta_t,0,1);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: ";
	tmp2 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 
	*/
	cout << "----------" << endl;
	cout << "Importance sampling, analytical local energy, numerical quantum force: " << endl;
	expectation = QD.Importance_Sampling_Metropolis_Expectation_Values(M,delta_t,1,0);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: ";
	tmp2 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 
	
	
	cout << "----------" << endl;
	cout << "Importance sampling, analytical local energy, analytical quantum force: " << endl;
	expectation = QD.Importance_Sampling_Metropolis_Expectation_Values(M,delta_t,1,1);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: ";
	tmp2 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 
	

	cout << "----------" << endl;
	cout << "Brute force, numerical local energy: " << endl;
	expectation = QD.Brute_Force_Metropolis_Expectation_Values(M,0);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: ";
	tmp2 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 

	
	cout << "----------" << endl;
	cout << "Brute force, analytical local energy: " << endl;
	expectation = QD.Brute_Force_Metropolis_Expectation_Values(M,1);
	cout << "Expectation_values:" << endl;
	cout << expectation << endl;
	cout << "variance: ";
	tmp2 = expectation(0);
	cout << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 
	
}