#ifndef project2lib_h
#define project2lib_h
#include <iostream>
#include <armadillo>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>
using namespace std;
using namespace arma;

class QuantumDots{
	/*
	The system.
	*/
	private:
		int N; //Number of electrons in trap.
	public:

		//Constructor
		QuantumDots(int a);

		//Print functions 
		void print_N_to_terminal(void);

};

#endif