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
	double result=0;

	
	mat phi_of_r_left =  zeros(N/2,N/2) ;
	mat phi_of_r_right =  zeros(N/2,N/2);
	for (int i = 0; i<=N-2; i += 2){
		for (int j=0; j<=N-2; j+=2){
			phi_of_r_left(j/2,i/2) = phi(i,r.col(j));
			phi_of_r_right(j/2,i/2) = phi(i+1,r.col(j+1));
		}
	}

	result = det(phi_of_r_left)*det(phi_of_r_right);
	double exponation_argument = 0;
	if (Jastrow_factor==1){
		mat r_i, r_j;
		double r_ij = 0;
		for (int j = 0; j<N; j++){
			for (int i = 0; i<j ; i++){
				r_i = r.col(i);
				r_j = r.col(j);
				r_ij = norm(r_i-r_j);
				exponation_argument += a(i,j)*r_ij/(1+beta*r_ij);
			}
		}
	}
	result *= exp(exponation_argument);
	


	//Analytical wavefunction 2 particle without repulsion
	/* 
	double r0 = norm(r.col(0));
	double r1 = norm(r.col(1));
	result = exp(-omega*alpha*(r0*r0 + r1*r1)/2.0);
	*/
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
	double tmp_wf = this->call(r);
	return pow(tmp_wf,2);
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
	double alpha_omega = alpha*omega;

	double fx = Hermite_polynomial(sqrt(alpha_omega)*x,nx(i));
	double fy = Hermite_polynomial(sqrt(alpha_omega)*y,ny(i));

	return fx*fy*exp(-alpha_omega*(x*x + y*y)/2);	
}

mat Trial_Wavefunction::nabla_phi(double l, int k, mat r_i){
	/*
	Function that returns nabla phi_k(r_i)/e^(-0.5*l2*r2) according to table 2.1.2.
	Input:
		- double l 			- The argument factor 
		- int k 	 		- funtion index
		- mat r_i 			- (2,1)-matrix, the position to evaluate the function. 
	Output:
		- mat result 		- A (2,1)-matrix with nabla phi_k(r_i)/e^(-0.5*l2*r2)
	*/
	mat result = zeros(2,1);

	double x = r_i(0,0), y = r_i(1,0);
	double x2 = x*x, y2 = y*y;
	double l2 = l*l;

	if (k==0){
		result(0,0) = -l2*x;
		result(1,0) = -l2*y;
	}
	else if (k==2){
		result(0,0) = -2*l*(l*x -1)*(l*x+1);
		result(1,0) = -2*l*l2*x*y;
	}
	else if(k==4){
		result(0,0) = -2*l*l2*x*y;
		result(1,0) = -2*l*(l*y - 1)*(l*y + 1);
	}
	else if (k==6){
		result(0,0) = -2*l2*x*(2*l2*x2 - 5);
		result(1,0) = -2*l2*y*(2*l2*x2 - 1);
	}
	else if (k==8){
		result(0,0) = -4*l2*y*(l*x-1)*(l*x+1);
		result(1,0) = -4*l2*x*(l*y-1)*(l*y+1);
	}
	else if (k==10){
		result(0,0) = -2*l2*x*(2*l2*y2 - 1);
		result(1,0) = -2*l2*y*(2*l2*y2 - 5);
	}
	return result;
}

double Trial_Wavefunction::nabla2_phi(double l,int k, mat r_i){
	/*
	Function that returns nabla phi^2_k(r_i) according to table 2.1.2.
	Input:
		- double l 			- The argument factor
		- int k 	 		- funtion index
		- mat r_i 			- (2,1)-matrix, the position to evaluate the function. 
	Output:
		- double result 	- nabla^2 phi_k(r_i)
	*/
	double result = 0;

	double x = r_i(0,0), y = r_i(1,0);
	double x2 = x*x, y2 = y*y;
	double l2 = l*l;
	double r2 = pow(norm(r_i),2);

	if (k==0){
		result = l2*(l2*r2-2);
	}
	else if (k==2){
		result = 2*l*l2*x*(l2*r2-4);
	}
	else if(k==4){
		result = 2*l*l2*y*(l2*r2-4);
	}
	else if (k==6){
		result = 2*l2*(l2*r2-6)*(2*l2*x2-1);
	}
	else if (k==8){
		result = 4*l2*l2*x*y*(l2*r2-6); 
	}
	else if (k==10){
		result = 2*l2*(l2*r2 - 6)*(2*l2*y2-1);
	}
	return result;
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
void QuantumDots::Set_Wavefunction(Trial_Wavefunction wf){
	/*
	Sets (updates) the wavefunction of the system.
	Input:
		- Trial_Wavefunction wf
	*/
	Wave_function = wf;
}
void QuantumDots::update_omega(double w){
	omega = w;
}



double QuantumDots::local_energy_function(mat r, int analytical_local_energy){
	/*
	Function that returns the local energy in a point of a wavefunction wf
	according to a hamiltonian operator H in a point r. 
	Input: 
		- double r 					- The point in which the local energy will be evaluated.
		- int R (default =0)		- Includes the electron-electron repulsion part of the Hamiltonian (1 if included, 0 else)
		- int analytical 		 	- To use or not to use analytical expressions
	Output:
		- double result 			- The local energy of the function. 
	*/

	//Initialization
	double result=0;
	double V = Potential(r);
	
	if (analytical_local_energy == 0){
		result = V - 0.5/Wave_function.call(r) *numerical_sum_laplacians(Wave_function,r);
	}
	else if(analytical_local_energy==1){
		//Laplacian term
		double scaled_laplacians = 0;
		for (int i=0;i<number_of_particles;i++){
			scaled_laplacians += LSP(Wave_function,r,i);
		}
		result = V - 0.5*scaled_laplacians;
	}
	else{
		cout << "-----------" << endl;
		cout << "Bad usage:" << endl;
		cout << "The 'analytical_local_energy' variable must be set to 1 or 0" << endl;
		cout << "Stopping program." << endl;
		cout << "-----------" << endl;
		exit(1);
	}
	
	return result;
}

//Help functions for the local energy
double QuantumDots::Potential (mat r){
	/*
	Function that evaluates the potential of a in a point r
	with the harmonic oscillator potential wit oscillator frequency omega (parameter of the class)
	with a electron repulsion part as well:
	Input:
		- mat r 						- Position in which the hamiltonian is to be evaluated
	Output:		
		-double V 						-The value of the potential in the point r.
	*/
	//Data

	double V = 0; 
	for (int i = 0; i<number_of_particles;i++){
		V += 0.5 * omega*omega*pow(norm(r.col(i)),2);
	}

	if (repulsion == 1){
		double r_ij = 0;
		for (int j = 0; j<number_of_particles; j++){
			for (int i = 0; i<j ; i++){
				r_ij = norm(r.col(i)-r.col(j));
				V += 1/r_ij;
			}
		} 
	}
	//cout << "Potiential: " << V << endl;
	return V;
}
double QuantumDots::Harmonic_potential(mat r){
	/*
	Function that returns just the harmonic part of the potential
	Input:
		- mat r 			- The position to be evaluated
	Output: 
		- double result 	- The harmonic potential at this point
	*/
	double result = 0; 
	for (int i = 0; i<number_of_particles;i++){
		result += 0.5 * omega*omega*pow(norm(r.col(i)),2);
	}
	return result;
}
double QuantumDots::Repulsive_potential(mat r){
	/*
	Function that returns just the repulsive potential 
	Input:
		- mat r 			- The position to be evaluated
	Output: 
		- double result 	- THe repulsive potential at the point.
	*/
	double result = 0;
	double r_ij = 0;
	if (repulsion==1){
		for (int j = 0; j<number_of_particles; j++){
			for (int i = 0; i<j ; i++){
				r_ij = norm(r.col(i)-r.col(j));
				result += 1/r_ij;
			}
		} 
	}

	//cout << "Potiential: " << result << endl;
	return result;
}

double QuantumDots::Analytical_kinetic_energy(mat r){
	/*
	Function that returns the local kinetic energy of the wavefunction.
	input: 
		- mat r 			- The position to be evaluated
	output:
		- double T 			- The local kinetic energy. 
	*/
	double T = 0;
	for (int i=0;i<number_of_particles;i++){
		T += LSP(Wave_function,r,i);
	}
	T *= -0.5;
	return T;
}


double QuantumDots::LSP(Trial_Wavefunction wf, mat r, int i){
	/*
	LSP - Laplacian sum part
	Function finding the analytical expressions NSS and N2SS (See theory) and returns the sum of these.
	Input:
		- Trial_Wavefunction wf 	- The wavefunction
		- mat r 					- The position in which to evaluate
		- int i 					- The particle number
	Output:
		- double result 			- N2SS + N2JJ + 2*NSS*NJJ 
	*/
	

	//Initialization	
	int N = wf.number_of_particles;
	double result=0;
	double r2 = pow(norm(r.col(i)),2);
	mat S_i;
	
	//Constructing Si
	if (i%2==0){
		S_i =  zeros(N/2,N/2) ;
		for (int it = 0; it<=N-2; it += 2){
			for (int j=0; j<=N-2; j+=2){
				S_i(j/2,it/2) = wf.phi(it,r.col(j));
			}
		}
	}
	else {
		S_i =  zeros(N/2,N/2);
		for (int it = 0; it<=N-2; it += 2){
			for (int j=0; j<=N-2; j+=2){
				S_i(j/2,it/2) = wf.phi(it+1,r.col(j+1));
			}
		}
	}
	mat S_i_inverse = inv(S_i);
	//S_i matrix is as wanted and so is the the inverse

	double N2SS=0;
	mat NSS = zeros(2,1);
	double l = sqrt(wf.omega*wf.alpha);
	for (int k=0;k<N/2;k++){
		NSS += S_i_inverse(k,i/2) * wf.nabla_phi(l,2*k,r.col(i)); 
		N2SS += S_i_inverse(k,i/2) * wf.nabla2_phi(l,2*k,r.col(i));
	}
	double exp_factor = exp(-0.5*l*l*r2);
	NSS *= exp_factor;
	N2SS *= exp_factor;

	mat NJJ = zeros(2,1);
	double N2JJ=0;
	double r_ik = 0, a_ik = 0, beta=0, a_ik_r_ik=0, denom=0;
	beta = wf.beta;
	for (int k=0;k<N;k++){
		if (k!=i){
			r_ik = norm(r.col(i)-r.col(k));
			a_ik = wf.a(i,k);
			a_ik_r_ik = a_ik/r_ik;
			denom =(1+beta*r_ik);

			NJJ += a_ik_r_ik*(r.col(i)-r.col(k))/(denom*denom);
			N2JJ += a_ik_r_ik*(1-beta*r_ik)/(denom*denom*denom);
		}
	}
	N2JJ += pow(norm(NJJ),2);

	result = N2SS;
	if (wf.Jastrow_factor==1){
		result += N2JJ + 2*dot(NJJ,NSS);
	}
	
	return result;
}

double QuantumDots::numerical_sum_laplacians(Trial_Wavefunction wf,mat r, double h ){
	/* ,
	Numerical function that returns the sum of det laplacians 
	in the position r = (v0 ... v_{N-1}) for each particle i acting on 
	a wavefunction object of type Trial_Wavefunction. 
	Input:
		- Trial_Wavefunction wf 		- Wavefunction to be evaluated
		- mat r 						- Point in which 
		- double h (default: 1e-4) 		- step length.
	Output: 
		-double result 				-The sum of the laplacians: sum_i (nabla_i^2) function
	*/
	int spatial_dimension = 2;
	int number_of_particles = wf.number_of_particles;
	double result = 0; //To store result


	mat r_jplus, r_jminus; //To store positions used in calc.
	double dfdj2_unscaled;  //To store double derivatives
	mat delta_r = zeros(spatial_dimension,number_of_particles);

	for (int i = 0; i<number_of_particles;i++){
		for (int j=0; j<spatial_dimension; j++){
			//dfdj2
			delta_r(j,i) = h;
			r_jminus = r - delta_r;
			r_jplus = r + delta_r;
			dfdj2_unscaled = wf.call(r_jplus) - 2*wf.call(r) + wf.call(r_jminus);
			delta_r(j,i) = 0;	
			result += dfdj2_unscaled;
		}
	}
	result /= (h*h);
	return result; 
}	
double QuantumDots::Average_distance(mat r){
	/*
	Function that returns the average distance between column vectors in a matrix. 
	Input:
		- mat r  	- The position matrix
	Output:
		- double Average_distance.
	*/
	double AD = 0;
	int tmp_counter = 0;
	for (int i = 0; i< number_of_particles; i++){
		for (int j=0; j<i; j++){
			AD += norm(r.col(j)-r.col(i));
			tmp_counter += 1;
		}
	}
	AD /= (double)tmp_counter;
	return AD;
}


//Monte Carlo Simulation functions
vec QuantumDots::Brute_Force_Metropolis_Expectation_Values(int M, int analytical_local_energy){
	/*
	Function that uses the Monte Carlo Metropolis algorithm to 
	find the the expectation value of the local energy and the local energy squared.
	Input:
		- int M 				 					- Number of Monte Carlo simulations
		- int analytical_local_energy (default 0)	- To use or not to use analytical expressions
	Output: 
		- vec expectation_values (2)- Expectation values of the local energy and the loc.energy squared
	*/
	//Extract relevant data 
	int number_of_particles = Wave_function.number_of_particles;

	//Initial position
	mat r = randn<mat>(2,number_of_particles);

	//Initialization
	int counter;
	int total_counter;
	double ratio;

	//Initial delta_r, will be modified to achieve 0.5 acceptance rate
	double delta_r;
	double acceptance_rate=2;


	//Find delta_r which gives around 0.5 acceptance rate: limits and number of test can be discussed. But 1000 works fine for now
	//Please note that this is actually already 25000 iterations, too much? 
	//Reduced to min(1/10 of MCS iterations,1000) 25.11.14 and verified that this gave almost the same acceptance rate and results.

	double wavefunction_squared = Wave_function.call_squared(r);
	for (double delta_ri = 0.2; delta_ri < 5; delta_ri += 0.1){
		counter = 0;
		total_counter = 0;
		for (int iter = 0; iter < min(M/10,1000); iter ++){
			//Choose particle to be moved
			int i = rand() % number_of_particles; 
			//Update step
			vec delta_vec_r = delta_ri * randn<vec>(2);
			mat r_p = r;
			r_p.col(i) = r_p.col(i) + delta_vec_r;
			double w = Wave_function.call_squared(r_p)/wavefunction_squared;
			double s = randu();
			if (w >= s){
				r = r_p;
				wavefunction_squared = Wave_function.call_squared(r);
				counter += 1;	
			}
			total_counter += 1;
		}
		double acceptance_rate_test = counter/(double)total_counter;
		if (abs(acceptance_rate_test-0.5) < abs(acceptance_rate-0.5)){
			acceptance_rate = acceptance_rate_test;
			delta_r = delta_ri;
		}
	}

	//Re-initialization before the real MCS
	counter = 0;
	total_counter = 0;
	r = randn<mat>(2,number_of_particles);
	wavefunction_squared = Wave_function.call_squared(r);
	double local_energy = local_energy_function(r,analytical_local_energy);
	double cumulative_local_energy = 0; 
	double cumulative_local_energy_squared = 0;
	
	
	while(total_counter < M){
		//Choose particle to be moved:
		int i = rand() % number_of_particles; 

		//Update step
		vec delta_vec_r = delta_r * randn<vec>(2);

		
		//Only coordenate move, x or y
		//int xory = rand() % 2;
		//delta_vec_r(xory,0) = 0;
		
	
		mat r_p = r;
		r_p.col(i) = r_p.col(i) + delta_vec_r;
		double w = Wave_function.call_squared(r_p)/wavefunction_squared;
		double s = randu();
		if (w >= s){
			r = r_p;
			wavefunction_squared = Wave_function.call_squared(r);
			counter += 1;	
			local_energy = local_energy_function(r,analytical_local_energy);

		}

		cumulative_local_energy += local_energy;
		cumulative_local_energy_squared += local_energy*local_energy;
		total_counter += 1;
	}

	//Uncomment if investigating ratio vs delta_r
	//ratio = counter/(double)(M-1);
	//cout << endl << "---- Acceptance rate: " << ratio;
	//cout << ", Step length: " << delta_r <<  "----" << endl;

	//Create matrix for storing expectation values
	vec expectation_values = zeros(2);
	expectation_values(0) = cumulative_local_energy/M;
	expectation_values(1) = cumulative_local_energy_squared/M;
	return expectation_values;
}	

vec QuantumDots::Importance_Sampling_Metropolis_Expectation_Values(int M, double delta_t, int analytical_local_energy, int analytical_quantum_force){
	/*
	Function that uses the Monte Carlo Metropolis algorithm to 
	find the the expectation value of the local energy and the local energy squared.
	Input:
		- int M 				 				- Number of Monte Carlo simulations
		- double delta_t 		 				- The time step
		- int analytical_local_energy 			- To use or not to use analytical expressions for the local energy
		- int analytical_quantum_force 			- To use or not to use analytical expressions for the quantum force
	Output: 
		- vec expectation_values (2)- Expectation values of the local energy and the loc.energy squared
	*/
	//Extract relevant data 
	int number_of_particles = Wave_function.number_of_particles;

	//Initial position
	mat r = randn<mat>(2,number_of_particles);

	//Initialization
	double local_energy = local_energy_function(r,analytical_local_energy);;
	double cumulative_local_energy = 0; 
	double cumulative_local_energy_squared = 0;
	int counter = 0;
	int total_counter = 0;
	double D = 0.5;
	double ratio;


	//Choose particle to be moved next:
	int i = rand() % number_of_particles; 

	mat quantum_force_old = quantum_force(r,analytical_quantum_force);
	double wavefunction_squared = Wave_function.call_squared(r);
	mat greensfunction_elements;
	double greensfunction;

	while(total_counter < M){
		
		//Update step

		vec delta_vec_r = D*quantum_force_old.col(i)*delta_t + sqrt(2*D*delta_t)*randn(2,1);


		/*
		//Only coordenate move, x or y
		int xory = rand() % 2;
		delta_vec_r(xory,0) = 0;
		*/

		mat r_p = r;
		r_p.col(i) = r_p.col(i) + delta_vec_r;

		mat quantum_force_new = quantum_force(r_p,analytical_quantum_force);
		double w = Wave_function.call_squared(r_p)/wavefunction_squared;
		double s = randu();


		//Greensfunction
		//My understanding 
		//greensfunction_elements = -0.25*(quantum_force_old + quantum_force_new)%
		//(2*(r_p.col(i)-r.col(i))/(D*delta_t) + quantum_force_new - quantum_force_old);

		//From the slides
		greensfunction_elements = 0.5*(quantum_force_old+quantum_force_new)%
		                          (D*delta_t*0.5*(quantum_force_old-quantum_force_new)-r_p + r);

		greensfunction = exp(sum(sum(greensfunction_elements))); 


		if (w*greensfunction >= s){
			r = r_p;
			wavefunction_squared = Wave_function.call_squared(r);
			counter += 1;	
			local_energy = local_energy_function(r,analytical_local_energy);
			quantum_force_old = quantum_force_new;
		}

		cumulative_local_energy += local_energy;
		cumulative_local_energy_squared += local_energy*local_energy;
		total_counter += 1;

		//Choose particle to be moved next:
		i = rand() % number_of_particles; 
	}
	
	//Uncomment if investigating ratio vs delta_r
	ratio = counter/(double)total_counter;
	cout << "Acceptance rate: " << ratio << endl;

	//Create matrix for storing expectation values
	vec expectation_values = zeros(2);
	expectation_values(0) = cumulative_local_energy/M;
	expectation_values(1) = cumulative_local_energy_squared/M;

	

	return expectation_values;
}	
vec QuantumDots::Metropolis_interesting_quantities(int M){
	/*
	Function that uses the Monte Carlo Metropolis algorithm to 
	find the the expectation value the energy, variance, average distance, potential_energy, harmonic potential energy,
	repulsive potential energy and kinetic_energy
	Input:
		- int M 				 					- Number of Monte Carlo simulations
	Output: 
		- vec expectation_values (7)- Expectation values in the order given above
	*/
	//Extract relevant data 
	int number_of_particles = Wave_function.number_of_particles;

	//Initial position
	mat r = randn<mat>(2,number_of_particles);

	//Initialization
	int counter;
	int total_counter;
	double ratio;

	//Initial delta_r, will be modified to achieve 0.5 acceptance rate
	double delta_r;
	double acceptance_rate=2;


	//Find delta_r which gives around 0.5 acceptance rate: limits and number of test can be discussed. But 1000 works fine for now
	//Please note that this is actually already 25000 iterations, too much? 
	//Reduced to min(1/10 of MCS iterations,1000) 25.11.14 and verified that this gave almost the same acceptance rate and results.
	double wavefunction_squared = Wave_function.call_squared(r);
	for (double delta_ri = 0.2; delta_ri < 5; delta_ri += 0.1){
		counter = 0;
		total_counter = 0;
		for (int iter = 0; iter < min(M/10,1000); iter ++){
			//Choose particle to be moved
			int i = rand() % number_of_particles; 
			//Update step
			vec delta_vec_r = delta_ri * randn<vec>(2);
			mat r_p = r;
			r_p.col(i) = r_p.col(i) + delta_vec_r;
			double w = Wave_function.call_squared(r_p)/wavefunction_squared;
			double s = randu();
			if (w >= s){
				r = r_p;
				wavefunction_squared = Wave_function.call_squared(r);
				counter += 1;	
			}
			total_counter += 1;
		}
		double acceptance_rate_test = counter/(double)total_counter;
		if (abs(acceptance_rate_test-0.5) < abs(acceptance_rate-0.5)){
			acceptance_rate = acceptance_rate_test;
			delta_r = delta_ri;
		}
	}

	//Re-initialization before the real MCS
	counter = 0;
	total_counter = 0;
	r = randn<mat>(2,number_of_particles);
	wavefunction_squared = Wave_function.call_squared(r);

	double harmonic_P = Harmonic_potential(r);
	double cumulative_HP = 0;

	double repulsive_P = Repulsive_potential(r);
	double cumulative_RP = 0;

	double potential_energy = harmonic_P + repulsive_P;
	double cumulative_PE = 0;

	double kinetic_energy = Analytical_kinetic_energy(r);
	double cumulative_KE = 0;

	double local_energy = potential_energy + kinetic_energy;
	double cumulative_local_energy = 0; 
	double cumulative_local_energy_squared = 0;

	double average_distance = Average_distance(r);
	double cumulative_AD = 0;
	
	
	while(total_counter < M){
		//Choose particle to be moved:
		int i = rand() % number_of_particles; 

		//Update step
		vec delta_vec_r = delta_r * randn<vec>(2);

		
		//Only coordenate move, x or y
		//int xory = rand() % 2;
		//delta_vec_r(xory,0) = 0;
		
	
		mat r_p = r;
		r_p.col(i) = r_p.col(i) + delta_vec_r;
		double w = Wave_function.call_squared(r_p)/wavefunction_squared;
		double s = randu();

		if (w >= s){
			r = r_p;
			wavefunction_squared = Wave_function.call_squared(r);
			counter += 1;
				harmonic_P = Harmonic_potential(r);
				repulsive_P = Repulsive_potential(r);
				potential_energy = harmonic_P + repulsive_P;
				kinetic_energy = Analytical_kinetic_energy(r);
				local_energy = potential_energy + kinetic_energy;
				average_distance = Average_distance(r);

		}
		cumulative_HP += harmonic_P;
		cumulative_RP += repulsive_P;
		cumulative_PE += potential_energy;
		cumulative_KE += kinetic_energy;
		cumulative_AD += average_distance;
		cumulative_local_energy += local_energy;
		cumulative_local_energy_squared += local_energy*local_energy;
		total_counter += 1;
		if (total_counter%100000==0){
			cout << "Progress: " << total_counter/(double)M*100  << " % " << endl;
		}
	}

	//Uncomment if investigating ratio vs delta_r
	//ratio = counter/(double)(M-1);
	//cout << endl << "---- Acceptance rate: " << ratio;
	//cout << ", Step length: " << delta_r <<  "----" << endl;

	//Create matrix for storing expectation values
	vec expectation_values = zeros(7);
	expectation_values(0) = cumulative_local_energy/M;
	expectation_values(1) = cumulative_local_energy_squared/M;
	expectation_values(2) = cumulative_AD/M;
	expectation_values(3) = cumulative_PE/M;
	expectation_values(4) = cumulative_HP/M;
	expectation_values(5) = cumulative_RP/M;
	expectation_values(6) = cumulative_KE/M;
	return expectation_values;
}	

//Importance sampling help function
mat QuantumDots::quantum_force(mat r, int analytical_quantum_force, double h){
	/*
	Returns the quantum force as described in the section on importance sampling. 
	Input:
		- mat r 						- The point in which the quantum force shall be returned
		- int i 						- The particle moved
		- int analytical_quantum_force 	- To use or not to use the analytical expressions
		- double h (Default = 1e-6)		- The step length in the nummerical evaluation
	Output:
		- vec result 					- The quantum force 1/psi del_i psi
	*/

	mat result = zeros(2,number_of_particles);

	if (analytical_quantum_force == 0){
		for (int i=0;i<number_of_particles;i++){
			mat rplusx = r;
			mat rminusx = r;
			rplusx(0,i) += h;
			rminusx(0,i) -= h;
			double dfdx = (Wave_function.call(rplusx)-Wave_function.call(rminusx))/(2*h);
			result(0,i) = dfdx;

			mat rplusy = r;
			mat rminusy = r;
			rplusy(1,i) += h;
			rminusy(1,i) -= h;
			double dfdy = (Wave_function.call(rplusy)-Wave_function.call(rminusy))/(2*h);
			result(1,i) = dfdy;
		}
		result *= 2/Wave_function.call(r);
	}
	else if (analytical_quantum_force == 1){
		for (int i=0;i<number_of_particles;i++){
			mat S_i;
			int N = number_of_particles;
		
			//Constructing Si
			if (i%2==0){
				S_i =  zeros(N/2,N/2) ;
				for (int it = 0; it<=N-2; it += 2){
					for (int j=0; j<=N-2; j+=2){
					S_i(j/2,it/2) = Wave_function.phi(it,r.col(j));
					}
				}
			}
			else {
				S_i =  zeros(N/2,N/2);
				for (int it = 0; it<=N-2; it += 2){
					for (int j=0; j<=N-2; j+=2){
						S_i(j/2,it/2) = Wave_function.phi(it+1,r.col(j+1));
					}
				}
			}
			mat S_i_inverse = inv(S_i);
			//S_i matrix is as wanted and so is the the inverse

			mat NSS = zeros(2,1);
			double l = sqrt(Wave_function.omega*Wave_function.alpha);
			for (int k=0;k<N/2;k++){
				NSS += S_i_inverse(k,i/2) * Wave_function.nabla_phi(l,2*k,r.col(i)); 
			}
			double r2 = pow(norm(r.col(i)),2);
			double exp_factor = exp(-0.5*l*l*r2);
			NSS *= exp_factor;

			mat NJJ = zeros(2,1);
			if (Wave_function.Jastrow_factor == 1){
				double r_ik = 0, a_ik = 0, beta=0, a_ik_r_ik=0, denom=0;
				beta = Wave_function.beta;
				for (int k=0;k<N;k++){
					if (k!=i){
						r_ik = norm(r.col(i)-r.col(k));
						a_ik = Wave_function.a(i,k);
						a_ik_r_ik = a_ik/r_ik;

						denom =(1+beta*r_ik);

						NJJ += a_ik_r_ik*(r.col(i)-r.col(k))/(denom*denom);
					}
				}
			}
			result.col(i) = 2*(NSS + NJJ);
		}
	}
	else{
		cout << "-----------" << endl;
		cout << "Bad usage:" << endl;	
		cout << "The 'analytical_quantum_force' variable must be set to 1 or 0" << endl;
		cout << "Stopping program." << endl;
		cout << "-----------" << endl;
		exit(1);
	}
	return result;
}



//Print functions 
void QuantumDots::print_numberofparticles_to_terminal(void){
	cout << "number_of_particles = " << number_of_particles << endl;
} 






















//********************Investigation class **************************//

//Constructors
Investigate::Investigate(double a0, double as, double am, double b0, double bs, double bm, QuantumDots sys){
	alpha_0 = a0; alpha_step = as; alpha_max = am;
	beta_0  = b0; beta_step  = bs; beta_max  = bm;

	alpha_dim = (alpha_max-alpha_0)/alpha_step;
	beta_dim = (beta_max-beta_0)/beta_step;

	system = sys;
}
Investigate::Investigate(QuantumDots sys){
	system = sys;
}


//Solve functions
void Investigate::find_minimum(int MCS,  int jastrow, int analytical_local_energy, int importance_sampling, int analytical_quantum_force, double delta_t){
	/*
	Function that finds the energies as functions of the parameters alpha and beta. 
	Stores the energies and the corresponding variances in  the matrices energies and variances.
	Input: 
		- int MCS				 		- Number of Monte Carlo simulations to be performed in finding each energy
		- int jastrow 					- With (jastrow=1) or without (jastrow=0) the jastrow factor.
		- int analytical_local energy 	- Numerical or analytical local energy expression
		- int importance sampling		- Using importance sampling (1) or brute force (0)
		- int analytical_quantum_force	- Using (1) the anlytical expression for the quantum force or not
		- double delta_t (default: 0.1)	- The time step in the importance sampling method. 

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


	//Find the dimensionality of the matrices
	int counter = 0;
	int total_counter = 0;
	double progress = 0;
	energies = zeros(beta_dim,alpha_dim);
	variances = zeros(beta_dim,alpha_dim);

	//Time mark
	clock_t start = clock();

	//progress bar
	
	cout << '|';
	for (int ai = 0; ai<alpha_dim; ai++){
		cout << '-';
	}
	cout << '|' << endl << '|' << flush;

	
	for (int ai = 0; ai<alpha_dim; ai++){
		counter += 1;
		#pragma omp parallel for num_threads(4)
		for (int bi = 0; bi<beta_dim; bi++){
			total_counter += 1;
			//alpha and beta
			double alpha = alpha_0 + ai*alpha_step;
			double beta = beta_0 + bi*beta_step;

			//Wavefunction
			Trial_Wavefunction wf (alpha,beta,omega,N,jastrow);
			QuantumDots copy (system.omega,system.number_of_particles,system.repulsion);
			copy.Set_Wavefunction(wf);

			vec expectation_values;
			//Energy and variance
			if (importance_sampling==0){
				expectation_values = copy.Brute_Force_Metropolis_Expectation_Values(MCS, analytical_local_energy);
			}
			else if (importance_sampling == 1){
				expectation_values = copy.Importance_Sampling_Metropolis_Expectation_Values(MCS,delta_t,analytical_local_energy,analytical_quantum_force);
			}	
			energies(bi,ai) = expectation_values(0);
			variances(bi,ai) = abs(expectation_values(1) - expectation_values(0)*expectation_values(0));
			
		}
		cout << '-' << flush;	

		/*
		//TIME
		progress = counter/(double)alpha_dim;
		//Time stamp
		clock_t now = clock();
		int factor = beta_dim;
		if (beta_dim>4){
			factor = 4;
		}

		cout << 100*progress << " %, " << "\t Time elapsed: " << (now-start)/(60*factor*(double)CLOCKS_PER_SEC) << " min ";
		double rate = (now-start)/((double)CLOCKS_PER_SEC*factor*progress);
		cout << "\t Time remaining: " << rate*(1-progress)/(double)60 << " min " << endl;
		*/

	}
	cout << '|' << endl;
	
}

void Investigate::find_parameters(int MCS, vec omegas_input, int jastrow){
	/*
	Function that finds the optimal parameteres for a range of omegas. 
	Uses the brute force approach with analytical local energy.
	Input:
		- MCS 						- Number of Monte Carlo simulations to be performed for each energy
		- vec omegas 				- Vector of omegas for which the parameters should be found
		- double resolution			- the resolution for the parameters
		- int jastrow 				- To have the jastrow factor or not
	Saves the optimal parameters in 
		-vec alpha_optimal
		-vec beta_optimal
		-vec energy_optimal
		-vec variance_optimal
	*/
	omegas = omegas_input;

	number_of_omegas = omegas.n_elem;
	alpha_optimal = zeros(number_of_omegas);
	beta_optimal = zeros(number_of_omegas);
	energy_optimal = zeros(number_of_omegas);
	variance_optimal=zeros(number_of_omegas);

	for (int index=0; index < number_of_omegas; index++){
		cout << "Progress: " << (index)/(double)number_of_omegas*100 << " %" << endl;
		system.update_omega(omegas(index));
		find_minimum(MCS,jastrow,1,0,0);
		uword alpha_index,beta_index;
		energy_optimal(index) = energies.min(beta_index,alpha_index);
		alpha_optimal(index) = alpha_0 + alpha_index*alpha_step;
		beta_optimal(index) = beta_0 + beta_index*beta_step;
		variance_optimal(index) = variances(beta_index,alpha_index);
	}
}



void Investigate::compare_analytical_numerical(int MCS, int jastrow){
	/*
	Function that investigates the difference between the analytical and numerical solutions as functions of alpha and beta. 
	Input: 
		- int MCS	 		- Number of Monte Carlo simulations to be performed in finding each energy
		- double delta_r 	- The step length to be used in the MC-simulation
		- int jastrow 		- With (jastrow=1) or without (jastrow=0) the jastrow factor.

	Saves <E_L>_analytical - <E_L>_numerical in matrix relative_energy_difference. 
	*/

	//Data
	int N = system.number_of_particles;
	double omega = system.omega;

	//Initialization
	double alpha, beta;
	vec expectation_values_analytical;
	vec expectation_values_numerical;


	//Find the dimensionality of the matrices
	int tot_counter = alpha_dim*beta_dim; int counter = 0;
	relative_energy_difference = zeros(beta_dim,alpha_dim);

	for (int ai = 0; ai<alpha_dim; ai++){
		#pragma omp parallel for num_threads(3) 
		for (int bi = 0; bi<beta_dim; bi++){
			counter += 1;
			//alpha and beta
			alpha = alpha_0 + ai*alpha_step;
			beta = beta_0 + bi*beta_step;

			//Wavefunction
			Trial_Wavefunction wf (alpha,beta,omega,N,jastrow);
			system.Set_Wavefunction(wf);

			//Energy and variance
			expectation_values_analytical = system.Brute_Force_Metropolis_Expectation_Values(MCS,1);
			expectation_values_numerical = system.Brute_Force_Metropolis_Expectation_Values(MCS,0);
			relative_energy_difference(bi,ai) = (expectation_values_analytical(0)-expectation_values_numerical(0))/expectation_values_analytical(0);
		}
	cout << counter/(double)tot_counter*100 << " %" << endl;
	}
}
void Investigate::analyze_importance_sampling(double alpha, double beta, int jastrow){
	/*
	Function that investigates the energies found with importance sampling against those found with 
	a brute force approach against both time interval dt and number of monte carlo simulations.
	Both with analytical expressions for the local energy and quantum force (for imp sampling).
	Input:
		- double alpha, beta 	- Trial wavefunction parameters
		- int jastrow 			- Whether to include the jastrow factor or not
	
	Saves the result in three matrices:
		- imp_dt 	 		- Time interval meshgrid
		- imp_MC 			- Monte Carlo simulations meshgrid
		- imp_energies		- Energies meshgrid for the importance sampling
		- brute_energies 	- Energies meshgrid for the brute force method

		Format:
					dt = 1e-6 	1e-5 	1e-4 	1e-3 	1e-2 	1e-1 	1 	10 	100
		MCS= 1e3 	E 			E 		... 									E
		MCS= 1e4 	E 			E 		... 									E
		MCS= 1e5 	...															...
		MCS= 1e6 	...															...
		MCS= 1e7 	E 			E 		... 									E
	*/
	//Initialization
	imp_dt = zeros(4,9);
	imp_MCS = zeros(4,9);
	imp_energies = zeros(4,9);
	brute_energies = zeros(4,9);

	Trial_Wavefunction wf (alpha,beta,system.omega,system.number_of_particles,jastrow);
	system.Set_Wavefunction(wf);

	//Loops over MCS first
	for (int MCS_index = 0; MCS_index < 4; MCS_index ++){
		int MCS = pow(10,MCS_index+3);

		#pragma omp parallel for 
		for (int dt_index = 0; dt_index < 9; dt_index ++){
			double dt = pow(10,dt_index-6);

			vec expect_brute = system.Brute_Force_Metropolis_Expectation_Values(MCS,1);
			vec expect_imp = system.Importance_Sampling_Metropolis_Expectation_Values(MCS,dt,1,1);

			brute_energies(MCS_index,dt_index) = expect_brute(0);
			imp_energies(MCS_index,dt_index) = expect_imp(0);
			imp_dt(MCS_index,dt_index) = dt;
			imp_MCS(MCS_index,dt_index) = MCS;
		}
		cout << "MCS = " << MCS << endl; 
	}

}

void Investigate::compare_times(double alpha, double beta, int jastrow){
	/*
	Program that compares the time used by the different algorithms to find the expectaiton value of the local energy.
	For different number of MCS (monte carlo simulations.) 
	Parallelization is not implemented since this was not demanded by the exercises.
	Input:
		- double alpha,beta 	- 	Parameters of the wavefunction used. 
		- int jastrow 			- Jastrow factor on or off

	Saves the time (in seconds) used in a matrix "times" of the following form 
	MCS-> 			10^3 	3*10^3	10^4 	3*10^4 	...		10^6
	Method 
	(BF,NLE)		t 		t 		... 					t
	(BF,ALE) 		t 		t 		...						t
	(IS,NLE,NQF)	...										...
	(IS,NLE,AQF)	...										...
	(IS,ALE,NQF)	...										...
	(IS,ALE,AQF)	t 		t 		... 					t

	And also the MCS meshgrid in "MCS_times"
	*/
	times = zeros(6,8);
	MCS_times = zeros(6,8);

	//Filling the MCS_times_matrix
	MCS_times(0,0) = pow(10,3); MCS_times(0,1) = 3*pow(10,3);
	MCS_times(0,2) = pow(10,4); MCS_times(0,3) = 3*pow(10,4);
	MCS_times(0,4) = pow(10,5); MCS_times(0,5) = 3*pow(10,5);
	MCS_times(0,6) = pow(10,6); MCS_times(0,7) = 3*pow(10,6);
	for (int i=1;i<6;i++){
		MCS_times.row(i) = MCS_times.row(0);
	}

	//Time variables:
	clock_t start, end;
	vec expect;

	Trial_Wavefunction wf (alpha,beta,system.omega,system.number_of_particles,jastrow);
	system.Set_Wavefunction(wf);

	for (int MCS_index = 0; MCS_index<8; MCS_index++){
		int MCS = MCS_times(0,MCS_index);
		cout << "MCS: " << MCS << endl;

		//(BF,NLE)
		start = clock();
		expect = system.Brute_Force_Metropolis_Expectation_Values(MCS,0);
		end = clock();
		times(0,MCS_index) = (end-start)/(double)CLOCKS_PER_SEC;

		//(BF,ALE)
		start = clock();
		expect = system.Brute_Force_Metropolis_Expectation_Values(MCS,1);
		end = clock();
		times(1,MCS_index) = (end-start)/(double)CLOCKS_PER_SEC;

		//(IS,NLE,NQF)
		start = clock();
		expect = system.Importance_Sampling_Metropolis_Expectation_Values(MCS,0.1,0,0);
		end = clock();
		times(2,MCS_index) = (end-start)/(double)CLOCKS_PER_SEC;

		//(IS,NLE,AQF)
		start = clock();
		expect = system.Importance_Sampling_Metropolis_Expectation_Values(MCS,0.1,0,1);
		end = clock();
		times(3,MCS_index) = (end-start)/(double)CLOCKS_PER_SEC;

		//(IS,ALE,NQF)
		start = clock();
		expect = system.Importance_Sampling_Metropolis_Expectation_Values(MCS,0.1,1,0);
		end = clock();
		times(4,MCS_index) = (end-start)/(double)CLOCKS_PER_SEC;

		//(IS,ALE,AQF)
		start = clock();
		expect = system.Importance_Sampling_Metropolis_Expectation_Values(MCS,0.1,1,1);
		end = clock();
		times(5,MCS_index) = (end-start)/(double)CLOCKS_PER_SEC;
	}
}
void Investigate::interesting_quantities(double alpha, double beta){
	/*
	Function that takes two trial function parameters alpha and beta and computes, 
	using the brute force metropolis algorithm with analytical expression for the local energy
	and  10^7 MCS, 
	several interesting quantities that are saved in the variables: 
	IQ_energy 				- The expectation value of the hamilton operator
	IQ_variance, 			- The variance of the local energy
	IQ_average_distance		- expectation value of the distance between each electron 
	IQ_potential_energy		- the expectation value of the potential
	IQ_harmonic_potential 	- The part of the potential due to harmonic oscillation
	IQ_repulsive_potential 	- The part of the potential du to repulsion
	IQ_kinetic_energy		- The expectation value of the kinetic energy
	
	Input: 	
		- double alpha, beta 	- Trial wavefunction parameters
	*/

	int MCS = pow(10,7);

	Trial_Wavefunction wf(alpha, beta,system.omega,system.number_of_particles,system.repulsion);
	system.Set_Wavefunction(wf);

	vec IQ = system.Metropolis_interesting_quantities(MCS);
	IQ_energy = IQ(0);
	IQ_variance = IQ(1)-IQ_energy*IQ_energy;
	IQ_average_distance = IQ(2);
	IQ_potential_energy = IQ(3);
	IQ_harmonic_potential = IQ(4);
	IQ_repulsive_potential = IQ(5);
	IQ_kinetic_energy = IQ(6);
}









//Print functions
void Investigate::print_energies_to_file(string filename){
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
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
	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			output << variances(bi,ai);
			if (bi != beta_dim-1){output << " ,";}
			}
		output << '\n';
		}
	output.close();
}

void Investigate::print_relative_energy_difference_to_file(string filename){
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			output << relative_energy_difference(bi,ai);
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
	for (int ai = 0; ai<alpha_dim; ai++){
		for (int bi = 0; bi<beta_dim; bi++){
			output << beta_0 + bi*beta_step;
			if (bi != beta_dim-1){output << " ,";}
			}
		output << '\n';
		}
	output.close();
}

void Investigate::print_optimals_to_file(string filename){
	/*
	Prints to file in the following manner
	omega0	alpha0 	beta0 	energy0 	variance0
	*/
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
	for (int index=0; index<number_of_omegas;index++){
		output << omegas(index) << ',';
		output << alpha_optimal(index)<<',';
		output << beta_optimal(index)<<',';
		output << energy_optimal(index)<<',';
		output << variance_optimal(index)<<',';
		if (index != number_of_omegas-1){
			output<< '\n';
		}
	}
	output.close();
}

void Investigate::print_importance_analysis_to_files(string filenamebase){
	/*
	Prints the results from the importance analysis to four csv-files:
	filenamebase_MCS_meshgrid, filenamebase_dt_meshgrid
	filenamebase_brute_energies, filenamebase_imp_energies,
	*/
	string 	fnb_MCS = filenamebase, 				fnb_dt = filenamebase, 				fnb_BE=filenamebase, 					fnb_IE=filenamebase;
			fnb_MCS.append("_MCS_meshgrid.csv"); 	fnb_dt.append("_dt_meshgrid.csv");	fnb_BE.append("_brute_energies.csv"); 	fnb_IE.append("_imp_energies.csv");

	const char* fnb_MCS_name = fnb_MCS.c_str();
	const char* fnb_dt_name = fnb_dt.c_str();
	const char* fnb_BE_name = fnb_BE.c_str();
	const char* fnb_IE_name = fnb_IE.c_str();

	ofstream MCS_out, dt_out, BE_out, IE_out;
	MCS_out.open(fnb_MCS_name);
	dt_out.open(fnb_dt_name);
	BE_out.open(fnb_BE_name);
	IE_out.open(fnb_IE_name);

	//Loops over MCS first
	for (int MCS_index = 0; MCS_index < 4; MCS_index ++){
		for (int dt_index = 0; dt_index < 9; dt_index ++){
			MCS_out << imp_MCS(MCS_index,dt_index) ;
			dt_out << imp_dt(MCS_index,dt_index) ;
			BE_out << brute_energies(MCS_index,dt_index) ;
			IE_out << imp_energies(MCS_index,dt_index) ;

			if (dt_index != 8){
				MCS_out << ',';
				dt_out << ',';
				BE_out << ',';
				IE_out << ',';
			}
		}
		if (MCS_index != 4){
				MCS_out << '\n';
				dt_out << '\n';
				BE_out << '\n';
				IE_out << '\n';
		}
	}
	MCS_out.close();
	dt_out.close();
	BE_out.close();
	IE_out.close();
}	


void Investigate::print_times_to_file(string filenamebase){
	/*
	Prints the result from the time analysis to two CSV-files
	filenamebase_times.csv and filenamebase_MCS.csv
	*/
	string 	fnb_MCS = filenamebase, 				fnb_times = filenamebase;			
			fnb_MCS.append("_MCS.csv"); 			fnb_times.append("_times.csv");	

	const char* fnb_MCS_name = fnb_MCS.c_str();
	const char* fnb_times_name = fnb_times.c_str();

	ofstream MCS_out, times_out;	
	MCS_out.open(fnb_MCS_name);
	times_out.open(fnb_times_name);

	//Loop over methods
	for (int method_index = 0; method_index<6; method_index++){
		for (int MCS_index = 0; MCS_index<8; MCS_index++){
			MCS_out << MCS_times(method_index,MCS_index);
			times_out << times(method_index,MCS_index);
			if (MCS_index != 7){
				MCS_out << ',';
				times_out << ',';
			}
		}
		if (method_index != 5){
				MCS_out << '\n';
				times_out << '\n';
		}
	}
	
	MCS_out.close();
	times_out.close();
}

void Investigate::print_interesting_quantities_to_file(string filename){
	/*
	Function that prints the interesting quantities to file. 
	Input:
		- string filename 	-  The name of the file to which the information is to be printed.
	*/
	ofstream output;
	const char* name = filename.c_str();
	output.open(name);
	output << "Interesting quantities of the wavefunction: " << endl;
	output << "----------"<< endl << endl;

	output << "Parameters: " << endl;
	output << "alpha: " << system.Wave_function.alpha << endl;
	output << "beta: " << system.Wave_function.beta << endl;
	output << "omega: " << system.omega << endl;
	output << "Number of particles " << system.number_of_particles << endl;
	output << "Number of MCS: 10^7" << endl;
	output << "Jastrow_factor: on" << endl;
	output << "----------" << endl << endl;

	output << "Expectation values: " << endl;
	output << "Energy:  " << IQ_energy<< endl;
	output << "Variance: " << IQ_variance<< endl;
	output << "Average distance: " << IQ_average_distance<< endl;
	output << "Potential energy: " << IQ_potential_energy<< endl;
	output << "Harmonic potential: " << IQ_harmonic_potential << endl;
	output << "Repulsive potential: " << IQ_repulsive_potential << endl;
	output << "Kinetic energy: " << IQ_kinetic_energy<< endl;

	output.close();
}






















//**********Functions needed for class and elsewhere**********//









//Test functions

mat twobytwo(void){
	mat A = zeros(2,2);
	A(0,0) = 0.84;
	A(1,0) = 0.39;
	A(0,1) = 0.78;
	A(1,1) = 0.79;
	return A;
}

mat twobysix(void){
	mat A = zeros(2,6);
	A(0,0) = 0.84;
	A(1,0) = 0.39;
	A(0,1) = 0.78;
	A(1,1) = 0.79;
	A(0,2) = 0.95;
	A(1,2) = 1.32;
	A(0,3) = 0.34;
	A(1,3) = 0.78;
	A(0,4) = 0.28;
	A(1,4) = 0.55;
	A(0,5) = 0.48;
	A(1,5) = 0.63;
	return A;
}
