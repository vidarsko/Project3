#include "project3lib.h"


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
	int total_counter = 0;

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
		total_counter += 1;
	}
	float ratio = counter/(float)total_counter;
	cout << ratio << endl;
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

	double result = exp(-omega*(r02 + r12)/2);
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
	a multi variable function f(r). 
	Input:
		- double (*f) 				- Function whose laplacian should be found
		- mat r 					- Point in which 
		- double h (default: 1e-6)	- step length.
	Output: 
		-double result 				-The sum of the laplacians: sum_i (nabla_i^2) function
	*/
	int spatial_dimension = r.n_rows;
	int number_of_particles = r.n_cols;
	double h2 = h*h;

	//Initialization
	double result = 0; //To store result
	mat r_jplus, r_jminus; //To store positions used in calc.
	double dfdj2;  //To store double derivatives
	mat delta_r = zeros(spatial_dimension,number_of_particles);

	for (int i = 0; i<number_of_particles;i++){
		for (int j=0; j<spatial_dimension; j++){
			//dfdj2
			delta_r(j,i) = h;
			r_jminus = r - delta_r;
			r_jplus = r + delta_r;
			dfdj2 = (f(r_jplus) - 2*f(r) + f(r_jminus))/h2;
			delta_r(j,i) = 0;	
			result += dfdj2;
		}
	}
	return result; 
}	

double laplacian_sum_ground_state(mat r){
	double r02 = pow(norm(r.col(0)),2);
	double r12 = pow(norm(r.col(1)),2);
	double omega = 1;
	double result;
	result = (omega*omega*(r02 + r12)-(4*omega))*two_particle_ground_state(r);
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
	int spatial_dimension = r.n_rows;
	int number_of_particles = r.n_cols;


	//The potential term:
	double V = 0; 
	for (int i = 0; i<number_of_particles;i++){
		V += 0.5 * omega*omega*pow(norm(r.col(i)),2);
	}
	double result = V*wf(r) - 0.5*sum_laplacians(wf,r);
	return result;
}





//Hermite polynomials
double Hermite_polynomial(double x, int degree){
	/*
	Function to evaluate the Hermite polynomial function.
	Input:
		- double x		- Point in which to evaulate the hermite polynomial
		- int degree	- The degree of the polynomial
	Output:
		-H(x) 			- The value of the polynomial at x
	*/
	double Hm1 = 1;
	double H = 2*x;
	double Hp1;
	for (int i = 0; i < degree;i++){
		Hp1 = 2*x*H - 2*(i+1)*Hm1;
		Hm1 = H;
		H = Hp1;
	}
	return Hm1;
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