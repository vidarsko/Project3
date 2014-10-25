#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	QuantumDots QD (1.0, 2, 1);
	Investigate inv(0.3, 0.1, 3.0,  0.001, 0.005, 0.1,   QD   );
	inv.solve(100000,2.5,1);
	inv.print_energies_to_file("test_energies.csv");
	inv.print_variances_to_file("test_variances.csv");
	inv.print_alpha_meshgrid_to_file("test_alpha.csv");
	inv.print_beta_meshgrid_to_file("test_beta.csv");
}