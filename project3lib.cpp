#include "project3lib.h"

double parameter = 1;

//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(double (*P)(mat), double (*g)(mat), int M, double delta_r, mat r){
	/*
	Function that uses the Monte Carlo Metropolis algorithm to 
	find the expectation value of some function g(r) with with 
	respect to a probabiltiy density P(r).
	Input:
		- double P 					- Funciton propto Prob density P
		- double g					- Function whose exp. value to be evaluated.
		- int M 					- Number of Monte Carlo simulations
		- double delta_r 			- Predefined step length
		- mat r 					- Initial position r
	Output: 
		- vec expectation_values (2)- Expectation values of g and g^2. 
	*/

	//Initialization
	double cumulative_function = 0; 		//Corresponds to local energy in the quantum example
	double cumulative_function_squared = 0;
	int counter = 0;

	while(counter < M){
		int i = rand() % r.n_cols; 
		vec delta_vec_r = delta_r * (randu<vec>(2)-0.5);
		mat r_p = r;
		r_p.col(i) = r_p.col(i) + delta_vec_r;
		vec tmp = randu<vec>(1); double s = tmp(0);
		double w = P(r_p)/P(r);
		if (w >= s){
			r = r_p;
			double gi = g(r);
			cumulative_function += gi;
			cumulative_function_squared += gi*gi;
			counter += 1;
		}
	}

	//Create matrix for storing expectation values
	vec expectation_values = zeros(2);
	expectation_values(0) = cumulative_function/M;
	expectation_values(1) = cumulative_function_squared/M;
	return expectation_values;
}	



//Test functions for the metropolis algorithm 
double Test_Probability_Density(mat r){
	/*
	Probability density where P(r)=e^(-radius) where radius = |mat r|:
	*/
	double radius = norm(r,2);
	return exp(-radius);
}
double Test_Evaluation_Function(mat r){
	/*
	Function to be evaluated in the probability density. 
	Here, the value is 
	*/
	return norm(r,2);
}
double two_particle_ground_state(mat r){
	/*
	The ground state function of the two particle system for the unperturbed 
	hamiltonian.
	Input:
		- mat r 						- The position to evaluate the function in.
		- double omega(default=1)		- The oscillator frequency
	Output: 
		- double result 				- The value of the wavefunction in that position
	*/
	double omega = 1;
	double r02 = pow(norm(r.col(0)),2);
	double r12 = pow(norm(r.col(1)),2);

	double result = exp(-parameter*omega*(r02 + r12)/2);
	return result; 
}



//Local energy evaluation function
double local_energy(double (*H)(double (*wf)(mat),mat),double (*wf)(mat),mat r){
	/*
	Function that returns the local energy in a point of a wavefunction wf
	according to a hamiltonian operator H in a point r. 
	Input: 
		- double (*wf) 		- The wavefunction
		- double (*H)		- The function evaluating the Hamiltonion of a wf in a point r
		- double r 			- The point in which the local energy will be evaluated.
	Output:
		- double le 		- The local energy of the function. 
	*/
	double le = 1/wf(r) * H(wf,r);
	//cout << le << endl;
	return le;
}

//Laplacian functions
double sum_laplacians(double (*f)(mat),mat r,double h=1e-4){
	/* ,
	Numerical function that returns the sum of det laplacians 
	in the position r = (v0 ... v_{N-1}) for each particle i acting on 
	a multi variable function "function"(r). 
	Input:
		- double (*f) 				- Function whose laplacian should be found
		- mat r 					- Point in which 
		- double h (default: 1e-6)	- step length.
	Output: 
		-double result 				-The sum of the laplacians: sum_i (nabla_i^2) function
	*/
	int spacial_dimension = r.n_rows;
	int number_of_particles = r.n_cols;
	double h2 = h*h;

	//Initialization
	double result = 0; //To store result
	mat r_xplus, r_xminus, r_yplus, r_yminus; //To store positions used in calc.
	double dfdx2, dfdy2;  //To store double derivatives
	mat delta_r = zeros(spacial_dimension,number_of_particles);

	for (int i = 0; i<number_of_particles;i++){

		//df/dx^2
		delta_r(0,i) = h;
		r_xminus = r - delta_r;
		r_xplus = r + delta_r;
		dfdx2 = (f(r_xplus) - 2*f(r) + f(r_xplus))/h2;
		delta_r(0,i) = 0;	

		//df/dy^2
		delta_r(1,i) = h;
		r_yminus = r - delta_r;
		r_yplus = r + delta_r;
		dfdy2 = (f(r_yplus) - 2*f(r) + f(r_yplus))/h2;
		delta_r(1,i) = 0;	

		result += dfdx2 + dfdy2;
	}
	return result; 
}	

double laplacian_sum_ground_state(mat r){
	double r02 = pow(norm(r.col(0)),2);
	double r12 = pow(norm(r.col(1)),2);
	double omega = 1;
	double result;
	result = ((parameter*omega)*(parameter*omega)*(r02 + r12)-(4*parameter*omega))*two_particle_ground_state(r);
	return result;
}


//Hamilton operator functions
double Unperturbed_Harmonic_Oscillator_Hamiltonian(double (*wf)(mat), mat r){
	/*
	Function that evaluates the hamiltonion of a wavefunction wf in a point r
	with an unperturbed harmonic oscillator potential:
	Input:
		- double (*wf)(mat) 		- Wavefunction
		- mat r 					- Position in which the hamiltonian is to be evaluated
		- double omega (default=1)	- The frequency of the oscillator
	Output:
		-double result 				-The value of the hamiltonian in the point r.
	*/
	double omega = 1;

	//Data
	int spacial_dimension = r.n_rows;
	int number_of_particles = r.n_cols;


	//The potential term:
	double V = 0; 
	for (int i = 0; i<number_of_particles;i++){
		V += 0.5 * omega*omega*pow(norm(r.col(i)),2);
	}
	cout << laplacian_sum_ground_state(r) - sum_laplacians(wf,r)<< endl;
	double result = V*wf(r) - 0.5*laplacian_sum_ground_state(r);
	//cout << (V*wf(r) - 0.5*laplacian_sum_ground_state(r))/wf(r) << endl;
	return result;

}














//**********************Quantum Dots class *******************//

//Constructor
QuantumDots::QuantumDots(int a){
	N = a;
}

//Print functions 
void QuantumDots::print_N_to_terminal(void){
	cout << "N = " << N << endl;
}












//******************Trial Wavefunction class ***************//

//Constructor
Trial_Wavefunction::Trial_Wavefunction(double a, double b){
	alpha = a;
	beta = b;
}