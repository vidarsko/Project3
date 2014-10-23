#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

double wf2(mat r){
		return pow(two_particle_ground_state(r),2);
}

	double le(mat r){
		return local_energy(Unperturbed_Harmonic_Oscillator_Hamiltonian,two_particle_ground_state,r);	
}

int main(){
	//Number of Monte Carlo simulations
	int M = 10;
	//Step length
	double delta_r = 2;
	//INitial position
	mat r = zeros(2,2);

	
	vec expect = Metropolis_Expectation_Values(wf2,le,M,delta_r,r);
	cout << "M = " << M << endl;
	cout << expect << endl;
	cout << "Variance: " << sqrt(expect(1)-expect(0)*expect(0))<< endl;


}