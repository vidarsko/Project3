#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	double alpha = 1, beta = 0, omega = 1;
	int number_of_particles = 12, spatial_dimension=2;
	Trial_Wavefunction wf (alpha,beta,omega,spatial_dimension,number_of_particles);
	QuantumDots QD(omega,number_of_particles);
	QD.Set_Wavefunction(wf);

	//Number of Monte Carlo simulations
	int M = 10000;
	//Step length
	double delta_r = 2.5; //Gives around 50% acceptance rate for N=2
	
	cout << "M = " << M << endl;
	cout << "alpha = " << alpha << endl;
	cout << "delta_r = " << delta_r << endl;
	vec expectation = Metropolis_Expectation_Values(QD,M,delta_r,1);
	cout << expectation;
	cout << "variance: " << abs(expectation(1)-expectation(0)*expectation(0)) << endl; 

}