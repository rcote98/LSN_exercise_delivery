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

int discrete_choice(double x){

	//return 1;

	x = x*6;

	if(x <= 1.){return 1;};
	if(x > 1  && x <= 2){return -1;};
	if(x > 2  && x <= 3){return 2;};
	if(x > 3  && x <= 4){return -2;};
	if(x > 4  && x <= 5){return 3;};
	if(x > 5){return -3;};

	return 1;
};

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
	// DISCRETE RANDOM WALK  ----------------------------------------

	unsigned int sims = 10000;        // Number of random walk simulations
	unsigned int steps = 100;       // Number of steps in each radnom walk	
	unsigned int a = 1;             // step size

	unsigned int M = sims * steps;   // Total number of throws
	unsigned int N = steps;           // Number of blocks
	
	double* r = new double[M];
	for(unsigned int i = 0; i < M; i++){	
		r[i] = rnd.Rannyu(); 
	}

	double ave[N];       // Averages vector
	double ave2[N];      // Squared Averages Vector
	double sum_prog[N];  //	Cumulative average
	double su2_prog[N];  // Cumulative square average
	double err_prog[N];  // Statistical uncertainty

	double pos[sims][3];
	for(unsigned int i = 0; i < sims; i++){
		pos[i][0] = 0;
		pos[i][1] = 0;
		pos[i][2] = 0;
	}



	int move = 0;
	double sum;
	unsigned int k;
	for(unsigned int i = 0; i<steps; i++){
		sum = 0;
		for(unsigned int j = 0; j<sims; j++){
			k = j + i*sims;

			move = discrete_choice(r[k]);
			if(move>0){pos[j][move-1] += a;}
			else{pos[j][abs(move)-1] -= a;}
			
			sum+= sqrt(pow(pos[j][0],2)+pow(pos[j][1],2)+pow(pos[j][2],2));
			
			
		}
		ave[i] = sum/sims;
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
	fout1.open("result2.2.discrete.csv", fstream::out);

	fout1 << "sim, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout1 << i+1 << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}
	fout1.close();

	// --------------------------------------------------------------
	// VECTOR RESET -------------------------------------------------

	for(unsigned int i = 0; i < N; i++){
		ave[i] = 0;
		ave2[i] = 0;
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		err_prog[i] = 0;
	}

	for(unsigned int i = 0; i < sims; i++){
		pos[i][0] = 0;
		pos[i][1] = 0;
		pos[i][2] = 0;
	}

	// ----------------------------------------------------------------
	// CONTINUOUS RANDOM WALK  ----------------------------------------

	double* rtheta = new double [M];
	double* rphi = new double[M];

	for(unsigned int i = 0; i < M; i++){	
		rtheta[i] = rnd.Angle()/2;
		rphi[i] = rnd.Angle();
	}
	
	for(unsigned int i = 0; i<steps; i++){
		sum = 0;
		for(unsigned int j = 0; j<sims; j++){
			k = j + i*sims;
			cout << k << ", " << rtheta[k] << ", " << rphi[k] << endl;

			pos[j][0] += a*sin(rtheta[k])*cos(rphi[k]);
			pos[j][1] += a*sin(rtheta[k])*sin(rphi[k]);
			pos[j][2] += a*cos(rtheta[k]);

			cout <<  pos[j][0] << ", " << pos[j][1]<< ", " <<  pos[j][2]  << endl;
			sum+= sqrt(pow(pos[j][0],2)+pow(pos[j][1],2)+pow(pos[j][2],2));
		}
		ave[i] = sum/sims;
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

	fout1.open("result2.2.continuous.csv", fstream::out);

	fout1 << "sim, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout1 << i+1 << ", " << sum_prog[i] << ", " << err_prog[i] 
		<< endl;
	}

	// -------------------------------------------------------------

	



	
	rnd.SaveSeed();
	return 0;
}



