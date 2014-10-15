#include "project3lib.h"




//**********************Quantum Dots class *******************//

//Constructor
QuantumDots::QuantumDots(int a){
	N = a;
}

//Print functions 
void QuantumDots::print_N_to_terminal(void){
	cout << "N = " << N << endl;
}