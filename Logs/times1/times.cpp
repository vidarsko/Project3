#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(20);

	//				omega 		N 		R
	QuantumDots QD (0.5, 		2, 		1);

	Investigate inv(QD);
	//					alpha 		beta 	jastrow					
	inv.compare_times(	0.77,		0.22,	1);
	inv.print_times_to_file("N2_alpha077_beta022");
}