#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	arma_rng::set_seed(20);

	// 				omega 	N 	R
	QuantumDots QD1(0.01, 2, 	0 );
	//
	Investigate inv1(QD1);
	//							alpha	beta
	inv1.interesting_quantities(1, 	0.07);
	inv1.print_interesting_quantities_to_file("N2norep_w001.txt");

	// 				omega 	N 	R
	QuantumDots QD2(0.1,  2, 	0 );
	//
	Investigate inv2(QD2);
	//							alpha	beta
	inv2.interesting_quantities(1, 	0.18);
	inv2.print_interesting_quantities_to_file("N2norep_w010.txt");

	// 					omega 	N 	R
	QuantumDots QD3(0.28, 2, 	0 );
	//
	Investigate inv3(QD3);
	//							alpha	beta
	inv3.interesting_quantities(1, 	0.26);
	inv3.print_interesting_quantities_to_file("N2norep_w028.txt");

	// 				omega 	N 	R
	QuantumDots QD4(0.5, 2, 	0 );
	//
	Investigate inv4(QD4);
	//							alpha	beta
	inv4.interesting_quantities(1, 	0.32);
	inv4.print_interesting_quantities_to_file("N2norep_w050.txt");

	// 				omega 	N 	R
	QuantumDots QD5(0.75, 2, 	0);
	//
	Investigate inv5(QD5);
	//							alpha	beta
	inv5.interesting_quantities(1, 	0.38);
	inv5.print_interesting_quantities_to_file("N2norep_w075.txt");

	// 				omega 	N 	R
	QuantumDots QD6(1, 2, 	0);
	//
	Investigate inv6(QD6);
	//							alpha	beta
	inv6.interesting_quantities(1, 	0.42);
	inv6.print_interesting_quantities_to_file("N2norep_w100.txt");


	// 				omega 	N 	R
	QuantumDots QD7(0.1, 6, 	0);
	//
	Investigate inv7(QD7);
	//							alpha	beta
	inv7.interesting_quantities(1, 	0.22);
	inv7.print_interesting_quantities_to_file("N6norep_w010.txt");

	// 				omega 	N 	R
	QuantumDots QD8(0.28,	6, 	0 );
	//
	Investigate inv8(QD8);
	//							alpha	beta
	inv8.interesting_quantities(1, 	0.33);
	inv8.print_interesting_quantities_to_file("N6norep_w028.txt");

	// 				omega 	N 	R
	QuantumDots QD9(0.5,	6, 	0);
	//
	Investigate inv9(QD9);
	//							alpha	beta
	inv9.interesting_quantities(1, 	0.38);
	inv9.print_interesting_quantities_to_file("N6norep_w050.txt");

	// 				omega 	N 	R
	QuantumDots QD10(0.75,	6, 	0 );
	//
	Investigate inv10(QD10);
	//							alpha	beta
	inv10.interesting_quantities(1, 	0.5);
	inv10.print_interesting_quantities_to_file("N6norep_w075.txt");

	// 				omega 	N 	R
	QuantumDots QD11(1,	6, 	0);
	//
	Investigate inv11(QD11);
	//							alpha	beta
	inv11.interesting_quantities(1, 	0.57);
	inv11.print_interesting_quantities_to_file("N6norep_w100.txt");

	// 				omega 	N 	R
	QuantumDots QD12(0.01,	6, 	0);
	//
	Investigate inv12(QD12);
	//							alpha	beta
	inv12.interesting_quantities(1, 	0.10);
	inv12.print_interesting_quantities_to_file("N6norep_w001.txt");

}