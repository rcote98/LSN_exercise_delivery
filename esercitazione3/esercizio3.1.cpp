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

// -------------------------------------------------------------------

double futureS(double t, double W, double S0, double mu, double sigma){
	return S0*exp((mu - 0.5*pow(sigma,2))*t  + sigma*W);
};

double stepS(double deltat, double Z, double S0, double mu, double sigma){
	return S0*exp((mu - 0.5*pow(sigma,2))*deltat  + sigma*Z*sqrt(deltat));;
}

// -------------------------------------------------------------------

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

	double t = 0;		 // starting time
	double S0 = 100;     // asset price at t = 0
	double T = 1;        // delivery time
	double K = 100;      // strike price
	double r = 0.1;      // risk-free interest rate
	double sigma = 0.25; // volatility

	// BLOCKING METHOD VARIABLES  -------------------------------

	unsigned int M = 100000;         // Total number of asset prices
	unsigned int N = 100;            // Number of blocks
	unsigned int L = int(M/N);       // Number of throws in each block

	double sum;
	unsigned int k;
	double *ave = new double[N];       // Averages vector
	double *ave2 = new double[N];      // Squared Averages Vector
	double *sum_prog = new double[N];  //	Cumulative average
	double *su2_prog = new double[N];  // Cumulative square average
	double *err_prog = new double[N];  // Statistical uncertainty

	fstream fout;

	// ----------------------------------------------------------
	// CALL DIRECT  -------------------------------------------

	double St, W, c_price;

	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			
			W = rnd.Gauss(0,T);

			St = futureS(T,W, S0, r, sigma);

			c_price = exp(-r*T)*max(0.,St-K);

			//c_price = call_option_price(t,St, sigma, r, K, T);

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
	
	fout.open("call_direct.csv", fstream::out);

	fout << "block, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << (i+1) << ", " << sum_prog[i] << ", " << err_prog[i] 
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
	

	unsigned int nstep = 100;
	double deltaT = T/nstep;

	double prevS, Z;

	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){

			prevS = S0;

			for(unsigned int k = 0; k<nstep; k++){
				Z = rnd.Gauss(0,1);
				prevS = stepS(deltaT, Z, prevS, r, sigma);
			}

			St = prevS;

			c_price = exp(-r*T)*max(0.,St-K);

			//c_price = call_option_price(t,St, sigma, r, K, T);

			
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
	
	fout.open("call_discrete.csv", fstream::out);
	fout << "block, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << (i+1) << ", " << sum_prog[i] << ", " << err_prog[i] 
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
	// PUT DIRECT  --------------------------------------------------

	double p_price;

	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			
			W = rnd.Gauss(0,T);

			St = futureS(T,W, S0, r, sigma);

			p_price = exp(-r*T)*max(0.,K-St);

			sum += p_price;
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
	
	fout.open("put_direct.csv", fstream::out);

	fout << "block, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << (i+1) << ", " << sum_prog[i] << ", " << err_prog[i] 
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
	// PUT DISCRETE  ------------------------------------------------

	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){

			prevS = S0;

			for(unsigned int k = 0; k<nstep; k++){
				Z = rnd.Gauss(0,1);
				prevS = stepS(deltaT, Z, prevS, r, sigma);
			}

			St = prevS;

			p_price = exp(-r*T)*max(0.,K-St);
			
			sum += p_price;
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
	
	fout.open("put_discrete.csv", fstream::out);

	fout << "block, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << (i+1) << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}
	fout.close();


	// --------------------------------------------------------------

	rnd.SaveSeed();
	return 0;

}
