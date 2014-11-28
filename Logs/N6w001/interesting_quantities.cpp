#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(20);

	// 					omega 	N 	R
	QuantumDots QD1(0.01, 6, 	1 );
	//
	Investigate inv1(QD1);
	//							alpha	beta
	inv1.interesting_quantities(0.61, 	0.10);
	inv1.print_interesting_quantities_to_file("N6w001.txt");


}