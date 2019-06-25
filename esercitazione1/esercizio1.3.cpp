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
	/*
		BUFFON'S PROBLEM SIMULATION

		Basically, the way the simulation works is by chosing a point
		p between 0 and d, and then sampling two random number x and y
		to choose a direction in 2D space without the need to use PI and
		then use that direction to project onto our original line.
		
		Then PI is calculated each block with the formula:
		PI = (2*Length)/(Probability*d) = (2*Throws*Length)/(Hits*d)

	*/	

	unsigned int M = 300000;         // Total number of throws
	unsigned int N = 100;            // Number of blocks
	unsigned int L = int(M/N);       // Number of throws in each block
	
	double* r = new double[3*M];	 // Random numbers to use
	for(unsigned int i = 0; i < 3*M; i++){	
		r[i] = rnd.Rannyu(); 
	}

	double ave[N];       // Averages vector
	double ave2[N];      // Squared Averages Vector
	double sum_prog[N];  //	Cumulative average
	double su2_prog[N];  // Cumulative square average
	double err_prog[N];  // Statistical uncertainty


	double d = 1.5;      // Distance between sticks
	double len = 1;      // Length of the sticks

	double p;            // random point between 0 and d
	double x, y;	     // random x, y points between [-1,1]
	double v;            // vector proyection ont original line
	
	double cpi;          // Calculated pi value
	double sum;          // Block sum for statistics

	unsigned int hits = 1; // Number of hits (starts in 1 to avoid division by zero)
	unsigned int k;        // Index used for the block structure
	
	for(unsigned int i = 0; i<N; i++){
		sum = 0;
		for(unsigned int j = 0; j < L; j++){
			
			k = j+i*L;
			
			p = r[k] * d;
			x = r[k +   M];
			y = r[k + 2*M];

			v = len * y/sqrt(pow(x,2) + pow(y,2));

			if(v>0){
				if(v > (d-p)) hits++;
			} else {
				if(v > p) hits++;
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
		fout1 << (i+1)*L << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}
	

	
	// ----------------------------------------------------------


	rnd.SaveSeed();
	return 0;

}
