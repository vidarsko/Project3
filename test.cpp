#include <iostream>
#include <armadillo>
#include "project3lib.h"
using namespace std;
using namespace arma;

int main(){
	QuantumDots test(2);
	int M = 10;
	double delta_r = 0.2;
	mat r = zeros(2,2);
	vec expect = Metropolis_Expectation_Values(Test_Probability_Density,Test_Evaluation_Function,M,delta_r,r);
	cout << "M = " << M << endl;
	cout << expect << endl;
	return 0;
}