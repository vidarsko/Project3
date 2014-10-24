#include "project3lib.h"






















//******************Trial Wavefunction class ***************//

//Constructor
Trial_Wavefunction::Trial_Wavefunction(double a, double b, double c, int S,int N){
	alpha = a;
	beta = b;
	omega = c;
	spatial_dimension = S;
	number_of_particles = N;
}
Trial_Wavefunction::Trial_Wavefunction(){};

double Trial_Wavefunction::call(mat r){
	/*
	Function that evaluates the multidimensional trial wavefunction at point r. 
	References to left and right matrices corresponds to the left and right 
	matrices in the evaluation of the modified slater determinant. 
	Jastrow factor not yet configured.
	Input:
		- mat r 		 	- The point at which the function is to be evaluated.
	Output: 
		- double result 	- The value of the wavefunction at this point. 
	*/
	int N = number_of_particles;

	mat phi_of_r_left =  zeros(N/2,N/2) ;
	mat phi_of_r_right =  zeros(N/2,N/2);
	for (int i = 0; i<=N-2; i += 2){
		for (int j=0; j<=N-2; j+=2){
			phi_of_r_left(j/2,i/2) = phi(i,r.col(j));
			phi_of_r_right(j/2,i/2) = phi(i+1,r.col(j+1));
		}
	}
	return det(phi_of_r_left)*det(phi_of_r_right);
}

double Trial_Wavefunction::call_squared(mat r){
	/*
	Function that evaluates the wavefunction squared at the point r. 
	I.e. the probability density.
	Input:
		- mat r 			- The point at which the function is to be evaluated.
	Output:
		- double result 	- The value of the squared wavefunction at that point.
	*/
	return pow(this->call(r),2);
}

//Help functions
double Trial_Wavefunction::Hermite_polynomial(double x, int degree){
	/*
	Function to evaluate the Hermite polynomial function based on the recursion formula.
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

double Trial_Wavefunction::nx(int i){
	/*
	nx(i) dependency, see theory.
	*/
	int result;
	if (i == 0 || i==1 || i==4 || i== 5 || i==10 || i== 11){result = 0;}
	else if (i==2 || i== 3 || i== 8 || i== 9){result = 1;}
	else if (i==6 || i== 7){result = 2;}
	return result;
}

double Trial_Wavefunction::ny(int i){
	/*
	ny(i) dependency, see theory.
	*/
	int result;
	if (i == 0 || i==1 || i==2 || i== 3 || i==6 || i== 7){result = 0;}
	else if (i==4 || i== 5 || i== 8 || i== 9){result = 1;}
	else if (i==10 || i== 11){result = 2;}
	return result;
}

double Trial_Wavefunction::phi(int i,mat r_i){
	/*
	Function that evaluates the trial function components phi_i(r_i) as described in theory.
	Input:
		- int i 		- Degree of the function component.
		- mat r_i 		- column matrix for which the component will be evaluated.
	Output: 
		- double result - The value of phi in r. 
	*/
	double x = r_i(0,0);
	double y = r_i(1,0);

	double fx = Hermite_polynomial(sqrt(alpha*omega)*x,nx(i))*exp(-alpha*omega*(x*x)/2);
	double fy = Hermite_polynomial(sqrt(alpha*omega)*y,ny(i))*exp(-alpha*omega*(y*y)/2);

	return fx*fy;	
}














































//**********************Quantum Dots class *******************//

//Constructor
QuantumDots::QuantumDots(int N){
	number_of_particles = N;
}

//Uncategorized
void QuantumDots::Set_Hamiltonian(double (*H) (Trial_Wavefunction wf, mat r)){
	/*
	Sets the Hamiltonian as a function of two paramters.
	Input:
		- double (*Hamiltonian) - Function
			Parameters of the hamiltonian:
				- Trial_Wavefunction Wave_function	- An wavefunction of type Trial_Wavefunction(object)
				- mat r 							- The point in which to evaluate the Hamiltonian
	*/
	Hamiltonian = H;
}

void QuantumDots::Set_Wavefunction(Trial_Wavefunction wf){
	/*
	Sets (updates) the wavefunction of the system.
	Input:
		- Trial_Wavefunction wf
	*/
	Wave_function = wf;
}


double QuantumDots::local_energy(mat r){
	/*
	Function that returns the local energy in a point of a wavefunction wf
	according to a hamiltonian operator H in a point r. 
	Input: 
		- double r 			- The point in which the local energy will be evaluated.
	Output:
		- double result 		- The local energy of the function. 
	*/
	double result = 1/Wave_function.call(r) * Hamiltonian(Wave_function,r);
	return result;
}

//Print functions 
void QuantumDots::print_numberofparticles_to_terminal(void){
	cout << "number_of_particles = " << number_of_particles << endl;
} 






















//**********Functions needed for class and elsewhere**********//

//Monte Carlo Simulation function
vec Metropolis_Expectation_Values(QuantumDots system, int M, double delta_r){
	/*
	Function that uses the Monte Carlo Metropolis algorithm to 
	find the the expectation value of the local energy and the local energy squared.
	Input:
		- QuantumDots system 		. The system (with wavefunction) to be evaluated
		- int M 					- Number of Monte Carlo simulations
		- double delta_r 			- Predefined step length
	Output: 
		- vec expectation_values (2)- Expectation values of g and g^2. 
	*/
	//Extract relevant data 
	int number_of_particles = system.Wave_function.number_of_particles;
	int spatial_dimension = system.Wave_function.spatial_dimension;

	//Initial position
	mat r = zeros(spatial_dimension,number_of_particles);

	//Initialization
	double local_energy = 0;
	double cumulative_local_energy = 0; 
	double cumulative_local_energy_squared = 0;
	int counter = 0;
	int total_counter = 0;

	while(counter < M){
		int i = rand() % number_of_particles; 
		vec delta_vec_r = delta_r * (randu<vec>(2)-0.5);
		mat r_p = r;
		r_p.col(i) = r_p.col(i) + delta_vec_r;
		vec tmp = randu<vec>(1); double s = tmp(0);
		double w = system.Wave_function.call_squared(r_p)/system.Wave_function.call_squared(r);
		cout << "w = " << w << " s = " << s << endl;
		if (w >= s){
			cout << "pling!" << endl;
			r = r_p;
			local_energy = system.local_energy(r);
			cumulative_local_energy += local_energy;
			cumulative_local_energy_squared += local_energy*local_energy;
			counter += 1;
		}
		total_counter += 1;
	}
	float ratio = counter/(float)total_counter;

	//Uncomment if investigating ratio vs delta_r
	//cout << ratio << endl;

	//Create matrix for storing expectation values
	vec expectation_values = zeros(2);
	expectation_values(0) = cumulative_local_energy/M;
	expectation_values(1) = cumulative_local_energy_squared/M;
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





//Laplacian functions
double sum_laplacians(Trial_Wavefunction wf,mat r,double h ){
	/* ,
	Numerical function that returns the sum of det laplacians 
	in the position r = (v0 ... v_{N-1}) for each particle i acting on 
	a wavefunction object of type Trial_Wavefunction. 
	Input:
		- Trial_Wavefunction wf 	- Wavefunction to be evaluated
		- mat r 					- Point in which 
		- double h (default: 1e-6)	- step length.
	Output: 
		-double result 				-The sum of the laplacians: sum_i (nabla_i^2) function
	*/
	int spatial_dimension = wf.spatial_dimension;
	int number_of_particles = wf.number_of_particles;
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
			dfdj2 = (wf.call(r_jplus) - 2*wf.call(r) + wf.call(r_jminus))/h2;
			delta_r(j,i) = 0;	
			result += dfdj2;
		}
	}
	return result; 
}	




//Hamilton operator functions
double Unperturbed_Harmonic_Oscillator_Hamiltonian(Trial_Wavefunction wf, mat r){
	/*
	Function that evaluates the hamiltonion of a wavefunction wf of type Trial_Wavefunction in a point r
	with an unperturbed harmonic oscillator potential wit oscillator frequency omega = 1:
	Input:
		- Trial_Wavefunction wf		- Wavefunction
		- mat r 					- Position in which the hamiltonian is to be evaluated
	Output:
		-double result 				-The value of the hamiltonian in the point r.
	*/
	double omega = 1;

	//Data
	int spatial_dimension = wf.spatial_dimension;
	int number_of_particles = wf.number_of_particles;


	//The potential term:
	double V = 0; 
	for (int i = 0; i<number_of_particles;i++){
		V += 0.5 * omega*omega*pow(norm(r.col(i)),2);
	}
	double result = V*wf.call(r) - 0.5*sum_laplacians(wf,r);
	return result;
}


