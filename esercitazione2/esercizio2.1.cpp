#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

// std function given the average, average squared and number of values n
double error(double av, double av2, int n){
	if (n==0) return 0;
	else return sqrt((av2 - av*av)/n);
};

// function to integrate
double funzione(double x){
	return M_PI*cos(M_PI*x/2)/2;
};

// taylor expansion of cos(x)
double taylor(double x){
  return 1.-pow(M_PI*x,2)/8.;
}

// distribution for importance sampling
double p_dist(double x){
  return (24./(24.- pow(M_PI, 2)))*taylor(x);
}

int main(int argc, char* argv[]){

	// RNG SETUP ----------------------------------------------

	Random rnd;

	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2;
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
	// CALCULATION WITHOUT IMPORTANCE SAMPLING ------------------

	unsigned int M = 100000;         // Total number of throws
	unsigned int N = 100;            // Number of blocks
	unsigned int L = int(M/N);       // Number of throws in each block

	double x[N];        // [0,1,2,..., N-1] * block size
	for(unsigned int i = 0; i < N; i++){
		x[i] = i+1;
	}

	double ave[N];       // Averages vector
	double ave2[N];      // Squared Averages Vector
	double sum_prog[N];  //	Cumulative average
	double su2_prog[N];  // Cumulative square average
	double err_prog[N];  // Statistical uncertainty

	double sum;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			sum += funzione(rnd.Rannyu());
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

	fstream fout1;
	fout1.open("result2.1.unif.csv", fstream::out);

	fout1 << "throws, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout1 << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}

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
	// CALCULATION WITH IMPORTANCE SAMPLING -------------------------
	
	double a, b;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){

			do{
				a = rnd.Rannyu();
				b = rnd.Rannyu();
			} while (b > p_dist(a));
			
			sum += funzione(a)/p_dist(a);

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
	
	fout1.close();
	fout1.open("result2.1.imp_sampling.csv", fstream::out);

	fout1 << "throws, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout1 << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}
	
	rnd.SaveSeed();
	return 0;

}
