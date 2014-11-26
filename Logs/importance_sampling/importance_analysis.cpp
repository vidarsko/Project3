#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(20);
	//				omega	N 	R
	QuantumDots QD (1.0, 	2, 	1);

	Investigate inv(QD);

	//								alpha 	beta 	jastrow
	inv.analyze_importance_sampling(0.72, 	0.24, 	1);

	// 										filenamebase
	inv.print_importance_analysis_to_files("imp_ana");
}