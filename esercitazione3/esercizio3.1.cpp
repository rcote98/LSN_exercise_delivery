#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

double error(double av, double av2, int n){
	if (n==0) return 0;
	else return sqrt((av2 - av*av)/n);
};


double N(double x){
	return (1. + erf(x/sqrt(2)))/2.;
};

double S_GBM_direct(double t, double S0, double mu, double sigma, double w){
	return S0*exp((mu - pow(sigma,2)/2.)*t + sigma*w);
}

double S_GBM_stepped(double t1, double t2, double S0, double mu, double sigma, double Z){
	return S0*exp((mu - pow(sigma,2)/2)*(t2-t1) + sigma*Z*sqrt(t2-t1));
}

double d1_func(double S, double r, double sigma, double T, double K, double t){
	return (log(S/K) + (r + pow(sigma, 2)*(T-t)/2))/(sigma*sqrt(T-t));
}

double d2_func(double d1, double sigma, double T, double t){
	return d1 - sigma*sqrt(T-t);
}

double call_price(double S, double r, double T, double K, double d1, double d2, double t){
	return S*N(d1) - K*exp(-r*(T-t))*N(d2);
}

double put_price(double S, double r, double T, double K, double d1, double d2, double t){
	return S*(N(d1) - 1) - K*exp(-r*(T-t))*(N(d2)-1);
}


int main(int argc, char* argv[]){

	// RNG SETUP ----------------------------------------------

	Random rnd;

	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	// ----------------------------------------------------------
	// Starting Variables

	double S0 = 100;     // asset price at t = 0
	double T = 1;        // delivery time
	double K = 100;      // strike price
	double r = 0.1;      // risk-free interest rate
	double sigma = 0.25; // volatility

	// BLOCKING METHOD VARIABLES  -------------------------------

	unsigned int M = 100000;         // Total number of asset prices
	unsigned int N = 100;            // Number of blocks
	unsigned int L = int(M/N);       // Number of throws in each block
	
	// Random numbers to use
	double* rand = new double[M*100]; // uniform between 0 and 1    	
	double* g = new double[M]; 

	for(unsigned int i = 0; i < M; i++){	
		g[i] = rnd.Gauss(0,T);
	}

	for(unsigned int i = 0; i < 100*M; i++){	
		rand[i] = rnd.Rannyu(); 
	}

	double x[N];        // [0,1,2,..., N-1] * block size
	for(unsigned int i = 0; i < N; i++){
		x[i] = i+1;
	}

	double sum;
	unsigned int k;
	double ave[N];       // Averages vector
	double ave2[N];      // Squared Averages Vector
	double sum_prog[N];  //	Cumulative average
	double su2_prog[N];  // Cumulative square average
	double err_prog[N];  // Statistical uncertainty

	fstream fout;

	// ----------------------------------------------------------
	// CALL DIRECT  -------------------------------------------


	double S, d1, d2, c_price;

	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			k = j+i*L;
			
			S = S_GBM_direct(T, S0, r, sigma, g[k]);
			d1 = d1_func(S, r, sigma, T, K, 0);
			d2 = d2_func(d1, sigma, T, 0);
			c_price = call_price(S, r, T, K, d1, d2, 0);
			sum += c_price;
		}
		ave[i] = sum/L;
		ave2[i] = pow(ave[i],2);
	}


	for(unsigned int i = 0; i<N; i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += ave2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}
	
	fout.open("result3.1.call_direct.csv", fstream::out);

	fout << "prize, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}
	fout.close();

	// --------------------------------------------------------------
	// VECTOR RESET -------------------------------------------------

	for(unsigned int i = 0; i < N; i++){
		ave[i] = 0;
		ave2[i] = 0;
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		err_prog[i] = 0;
	}

	// --------------------------------------------------------------
	// CALL DISCRETE  -----------------------------------------------


	unsigned int interval_number = 100;
	double deltaT = T/100; 
	double t1, t2;

	
	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			S = S0;
			t1 = -deltaT;
			t2 = 0;
			for (unsigned int p = 1; p <= interval_number; p++){
				k = (p-1) + j*interval_number+i*L;
				t1 += deltaT;
				t2 += deltaT;
				S = S_GBM_stepped(t1, t2, S, r, sigma, rand[k]);
			}

			d1 = d1_func(S, r, sigma, T, K, 0);
			d2 = d2_func(d1, sigma, T, 0);
			c_price = call_price(S, r, T, K, d1, d2, 0);
			
			sum += c_price;
		}
		ave[i] = sum/L;
		ave2[i] = pow(ave[i],2);
	}


	for(unsigned int i = 0; i<N; i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += ave2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}
	
	fout.open("result3.1.call_discrete.csv", fstream::out);
	fout << "prize, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}
	fout.close();




	rnd.SaveSeed();
	return 0;

}
