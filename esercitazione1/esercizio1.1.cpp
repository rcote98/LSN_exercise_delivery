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
	// VALORE MEDIO DI R ----------------------------------------

	unsigned int M = 100000;         // Total number of throws
	unsigned int N = 100;            // Number of blocks
	unsigned int L = int(M/N);       // Number of throws in each block
	
	double r[M];	    // Random numbers to use
	for(unsigned int i = 0; i < M; i++){	
		r[i] = rnd.Rannyu(); 
	}

	double x[N];        // [0,1,2,..., N-1]
	for(unsigned int i = 0; i < N; i++){
		x[i] = i;
	}

	double ave[N];       // Averages vector
	double ave2[N];      // Squared Averages Vector
	double sum_prog[N];  //	Cumulative average
	double su2_prog[N];  // Cumulative square average
	double err_prog[N];  // Statistical uncertainty

	double sum;
	unsigned int k;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			k = j+i*L;
			sum += r[k];
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
	
	for(unsigned int i = 0; i < N; i++){
		x[i] = x[i]*L;
	}

	fstream fout1;
	fout1.open("result1.1.r.csv", fstream::out);

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
	// VALORE MEDIO DI SIGMA

	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			k = j+i*L;
			sum += (r[k] - 0.5)*(r[k] - 0.5);
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

	fstream fout2;
	fout2.open("result1.1.sigma.csv", fstream::out);

	fout2 << "throws, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout2 << x[i] << ", " << sum_prog[i] << ", " 
		<< err_prog[i] << endl;
	}


	// -----------------------------------------------------------
	// TEST DI CHI QUADRO
	
	unsigned int n = 10000;     // Number of random numbers
	unsigned int T = 100;       // Number of intervals

	double random[n];           // Vector containing the random number
	double ni[T];               // Placeholder for count number

	double chi[100];            // Container for chi sqaured

	double min, max;
	for(unsigned int j = 0; j < 100; j++){

		for(unsigned int i = 0; i < n; i++){ni[i] = 0;}
		for(unsigned int i = 0; i < n; i++){random[i] = rnd.Rannyu();}

		

		for(unsigned int j= 0;j<T; j++){
			min = ((double) j)*1./T;
			max = ((double)j+1)*1/T;
			for(unsigned int i = 0; i < n; i++){
				if(random[i] > min && random[i] < max) ni[j] += 1;
			}
			
			sum = 0;
			for(unsigned int i = 0; i < T; i++){
				sum += pow((ni[j] - (n/T)), 2)/(n/T);
			}
			chi[j] = sum;
		} 
	}
	
	fstream fout3;
	fout3.open("result1.1.chi.csv", fstream::out);

	fout3 << "chij" << endl;
	for(unsigned int j = 0; j < 100; j++){
		fout3 << chi[j] << endl;
	}




	
	rnd.SaveSeed();
	return 0;

}
