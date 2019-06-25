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
	else return sqrt(pow(av-av2, 2)/n);
}

// dice choice out of a random number
int dice(double r){
	if(r < 1./6.) return 1;
	if(r > 1./6. && r < 2./6.) return 2;
	if(r > 2./6. && r < 3./6.) return 3;
	if(r > 3./6. && r < 4./6.) return 4;
	if(r > 4./6. && r < 5./6.) return 5;
	else return 6;
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
	

	unsigned int M = 10000;                    // Total number of throws
	unsigned int N[4] = {1,2,10,100};          // Number of throws
	
	double sum;	
	double series_sum[4][M];

	// GAUSSIAN -------------------------------------------------

	for(unsigned int i=0; i<4; i++){
		for(unsigned int j=0; j<M; j++){
			sum = 0;
			for(unsigned int k=0; k < N[i]; k++){
				sum += (double)dice(rnd.Rannyu())/(double)N[i];
			}
			series_sum[i][j] = sum;	
		}
	}

	fstream fout1;
	fout1.open("result1.2.dice.csv", fstream::out);
	fout1 << "N=1, N=2, N=10, N=100" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout1 << series_sum[0][i] << ", " << series_sum[1][i] << ", " 
		<< series_sum[2][i] << ", " << series_sum[3][i] << endl;
	}
	
	// VECTOR RESET -----------------------------------------------
	
	for(unsigned int i=0; i < 4; i++){
		for(unsigned int j=0; j < M; j++){
			series_sum[i][j] = 0;
		}
	}

	// -------------------------------------------------------------
	// EXPONENTIAL -------------------------------------------------

	for(unsigned int i=0; i<4; i++){
		for(unsigned int j=0; j<M; j++){
			sum = 0;
			for(unsigned int k=0; k < N[i]; k++){
				sum += rnd.Exponential(1)/N[i];
			}
			series_sum[i][j] = sum;	
		}
	}

	fstream fout2;
	fout2.open("result1.2.exp.csv", fstream::out);
	fout2 << "N=1, N=2, N=10, N=100" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout2 << series_sum[0][i] << ", " << series_sum[1][i] << ", " 
		<< series_sum[2][i] << ", " << series_sum[3][i] << endl;
	}
	
	// VECTOR RESET -----------------------------------------------
	
	for(unsigned int i=0; i < 4; i++){
		for(unsigned int j=0; j < M; j++){
			series_sum[i][j] = 0;
		}
	}

	// -------------------------------------------------------------
	// LORENTZIAN --------------------------------------------------	
	
	for(unsigned int i=0; i<4; i++){
		for(unsigned int j=0; j<M; j++){
			sum = 0;
			for(unsigned int k=0; k < N[i]; k++){
				sum += rnd.Lorentz(0,1)/N[i];
			}
			series_sum[i][j] = sum;	
		}
	}

	fstream fout3;
	fout3.open("result1.2.lorentz.csv", fstream::out);
	fout3 << "N=1, N=2, N=10, N=100" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout3 << series_sum[0][i] << ", " << series_sum[1][i] << ", " 
		<< series_sum[2][i] << ", " << series_sum[3][i] << endl;
	}

	rnd.SaveSeed();
	return 0;

}
