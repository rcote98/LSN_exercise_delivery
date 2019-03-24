#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

double error(double av, double av2, int n){
	if (n==0) return 0;
	else return sqrt(pow(av-av2, 2)/n);
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
	// pi = 2L/Pd
	unsigned int M = 300000;         // Total number of throws
	unsigned int N = 100;            // Number of blocks
	unsigned int L = int(M/N);       // Number of throws in each block
	
	double r[2*M];	    // Random numbers to use
	for(unsigned int i = 0; i < 2*M; i++){	
		r[i] = rnd.Rannyu(); 
	}

	double throws[N];
	for(unsigned int i = 0; i < N; i++){	
		throws[i] = (i+1)*L;
	}
	

	double ave[N];       // Averages vector
	double ave2[N];      // Squared Averages Vector
	double sum_prog[N];  //	Cumulative average
	double su2_prog[N];  // Cumulative square average
	double err_prog[N];  // Statistical uncertainty


	double d = 1.5;
	double len = 1;

	double p;
	double x;
	
	double cpi;
	double sum;
	unsigned int hits = 1;
	unsigned int k;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			k = j+i*L;
			
			p = d + r[k];
			x = -len + 2*len*r[k+M];

			if(x>0){
				if( fmod(p+x, d) < x) hits++;
			} else {
				if( d- (fmod(p+x, d) < x) hits++;
			};
			
			cpi = (2*len*k)/(hits*d);
			sum += cpi;;
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
	fout1.open("result1.3.pi.csv", fstream::out);

	fout1 << "throws, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout1 << throws[i] << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}
	

	
	// ----------------------------------------------------------


	rnd.SaveSeed();
	return 0;

}
