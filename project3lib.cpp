#include "project3lib.h"




//******************Trial Wavefunction class ***************//

//Constructors
Trial_Wavefunction::Trial_Wavefunction(double a, double b, double c,int N,int J){
	/*
	Trial Wavefunction constructor.
	Input:
		- double a, b 		- alpha, beta - variational parameters.
		- double c 			- omega - potential oscillattion frequency 	
		- int N 			- number of particles
		- int J 			- Jastrow factor
	*/
	alpha = a;
	beta = b;
	omega = c;
	number_of_particles = N;
	Jastrow_factor = J;
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
	double result;

	mat phi_of_r_left =  zeros(N/2,N/2) ;
	mat phi_of_r_right =  zeros(N/2,N/2);
	for (int i = 0; i<=N-2; i += 2){
		for (int j=0; j<=N-2; j+=2){
			phi_of_r_left(j/2,i/2) = phi(i,r.col(j));
			phi_of_r_right(j/2,i/2) = phi(i+1,r.col(j+1));
		}
	}

	result = det(phi_of_r_left)*det(phi_of_r_right);

	if (Jastrow_factor==1){
		mat r_i, r_j;
		double r_ij = 0;
		for (int j = 0; j<N; j++){
			for (int i = 0; i<j ; i++){
				r_i = r.col(i);
				r_j = r.col(j);
				r_ij = norm(r_i-r_j);
				result *= exp(a(i,j)*r_ij/(1+beta*r_ij));
			}
		}
	}
	return result;

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


double Trial_Wavefunction::a(int i, int j){
	/*
	Function that returns 1/3 if spin of i and j are parallell and 1 if they are anti-parallel
	(see theory on the Jastrow factor)
	Input:
		- int i, j 		 	-Particle indices
	Output: 
		- double 1 or 1/3 	-Factor in the Jastrow factor
	*/
	int test = i + j;
	double result;
	if (test%2==0){ 
		result = 1.0/3.0;
	}
	else {
		result = 1;
	}
	return result;
}















































//**********************Quantum Dots class *******************//

//Constructor
QuantumDots::QuantumDots(double w, int N, int R){
	omega = w;
	number_of_particles = N;
	repulsion = R;
}
QuantumDots::QuantumDots(){}


//Configuration member functions
double QuantumDots::Hamiltonian (Trial_Wavefunction wf, mat r){
	/*
	Function that evaluates the hamiltonion of a wavefunction wf of type Trial_Wavefunction in a point r
	with the harmonic oscillator potential wit oscillator frequency omega (parameter of the class)
	with a electron repulsion part as well:
	Input:
		- Trial_Wavefunction wf		- Wavefunction
		- mat r 					- Position in which the hamiltonian is to be evaluated
		- int R 					- Includes the electron-electron repulsion part of the Hamiltonian (1 if included, 0 else)
	Output:
		-double result 				-The value of the hamiltonian in the point r.
	*/

	//Data
	int number_of_particles = wf.number_of_particles;


	//The potential term:
	double V = 0; 
	for (int i = 0; i<number_of_particles;i++){
		V += 0.5 * omega*omega*pow(norm(r.col(i)),2);
	}
	double result = V*wf.call(r) - 0.5*sum_laplacians(wf,r);
	if (repulsion == 1){
		mat r_i, r_j;
		double r_ij = 0;
		for (int j = 0; j<number_of_particles; j++){
			for (int i = 0; i<j ; i++){
				r_i = r.col(i);
				r_j = r.col(j);
				r_ij = norm(r_i-r_j);
				result += 1/r_ij;
			}
		} 
	}
	return result;
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
		- double r 				 - The point in which the local energy will be evaluated.
		- int R (default =0)	 - Includes the electron-electron repulsion part of the Hamiltonian (1 if included, 0 else)
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






















//********************Investigation class **************************//

//Constructor
Investigate::Investigate(double a0, double as, double am, double b0, double bs, double bm, QuantumDots sys){
	alpha_0 = a0; alpha_step = as; alpha_max = am;
	beta_0  = b0; beta_step  = bs; beta_max  = bm;

	system = sys;
}

//Solve function
void Investigate::solve(int MCS, double delta_r, int jastrow){
	/*
	Function that finds the energies as functions of the parameters alpha and beta. 
	Stores the energies and the corresponding variances in  the matrices energies and variances.
	Input: 
		- int MCS	 		- Number of Monte Carlo simulations to be performed in finding each energy
		- double delta_r 	- The step length to be used in the MC-simulation
		- int jastrow 		- With (jastrow=1) or without (jastrow=0) the jastrow factor.

	Structure of energies and variances:
							alpha_0 	alpha_0+alpha_step  	alpha_0+2alpha_step 	... 	alpha_max-delta1
	beta_0 					E 			E 						E 						... 	E 
	beta_0+beta_step		E 			E 						E 						... 	E 
	beta_0+2beta_step		...
	...						...
	beta_max-delta2			E 			E 						E 						...		E
	*/

	//Data
	int N = system.number_of_particles;
	double omega = system.omega;

	//Initialization
	double alpha, beta;
	vec expectation_values = zeros(2);


	//Find the dimensionality of the matrices
	int alpha_dim = (alpha_max-alpha_0)/alpha_step;
	int beta_dim = (beta_max-beta_0)/beta_step;
	int tot_counter = alpha_dim*beta_dim; int counter = 0;
	energies = zeros(beta_dim,alpha_dim);
	variances = zeros(beta_dim,alpha_dim);

	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			counter += 1;
			//alpha and beta
			alpha = alpha_0 + ai*alpha_step;
			beta = alpha_0 + bi*beta_step;

			//Wavefunction
			Trial_Wavefunction wf (alpha,beta,omega,N,jastrow);
			system.Set_Wavefunction(wf);

			//Energy and variance
			expectation_values = Metropolis_Expectation_Values(system,MCS,delta_r);
			energies(bi,ai) = expectation_values(0);
			variances(bi,ai) = abs(expectation_values(1) - expectation_values(0)*expectation_values(0));

		}
	cout << counter/(double)tot_counter*100 << " %" << endl;
	}
}

void Investigate::print_energies_to_file(string filename){
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
	int alpha_dim = energies.n_cols;
	int beta_dim = energies.n_rows;
	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			output << energies(bi,ai);
			if (bi != beta_dim-1){output << " ,";}
			}
		output << '\n';
		}
	output.close();
}

void Investigate::print_variances_to_file(string filename){
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
	int alpha_dim = variances.n_cols;
	int beta_dim = variances.n_rows;
	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			output << variances(bi,ai);
			if (bi != beta_dim-1){output << " ,";}
			}
		output << '\n';
		}
	output.close();
}

void Investigate::print_alpha_meshgrid_to_file(string filename){
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
	int alpha_dim = energies.n_cols;
	int beta_dim = energies.n_rows;
	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			output << alpha_0 + ai*alpha_step;
			if (bi != beta_dim-1){output << " ,";}
			}
		output << '\n';
		}
	output.close();
}

void Investigate::print_beta_meshgrid_to_file(string filename){
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
	int alpha_dim = energies.n_cols;
	int beta_dim = energies.n_rows;
	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			output << beta_0 + bi*beta_step;
			if (bi != beta_dim-1){output << " ,";}
			}
		output << '\n';
		}
	output.close();
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
	int spatial_dimension = 2;

	//Initial position
	mat r = randu<mat>(spatial_dimension,number_of_particles);

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
		if (w >= s){
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
	int spatial_dimension = 2;
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




