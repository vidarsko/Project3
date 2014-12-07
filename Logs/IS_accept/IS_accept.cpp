#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(20);
	//				omega	N 	R
	QuantumDots QD (1.0, 	2, 	1);
	Trial_Wavefunction wf (0.72,0.24,1,2,1);
	QD.Set_Wavefunction(wf);

	for (int dt_index = -6; dt_index < 3; dt_index ++){
		double dt = pow(10,dt_index);
		cout << "dt = " << dt << endl;
		vec expect = QD.Importance_Sampling_Metropolis_Expectation_Values(pow(10,7), dt,1, 1);
	}
}