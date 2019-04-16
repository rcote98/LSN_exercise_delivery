/****************************************************************
*****************************************************************
		_/		_/	_/_/_/	_/			 Numerical Simulation Laboratory
	 _/_/  _/ _/			 _/				Physics Department
	_/	_/_/		_/		_/			 Universita' degli Studi di Milano
 _/		 _/				_/ _/				Prof. D.E. Galli
_/		_/	_/_/_/	_/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>			// srand, rand: to generate random number
#include <iostream>			// cin, cout: Standard Input/Output Streams Library
#include <fstream>			// Stream class to both read and write from/to files.
#include <cmath>            // Math stuff

#include "random.h"

using namespace std;

double error(double av, double av2, int n);
void data_blocking(unsigned int N, double* ave, string fname);

double* metropolis(double (wave)(double x, double y, double z), double* old_pos, double a, Random *rnd);

double base_state(double x, double y, double z);
double excited_state(double x, double y, double z);



int main(){
 
	// ----------------------------------------------------------------
	// RNG SETUP
	// ----------------------------------------------------------------

	Random *rnd = new Random();

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
				rnd->SetRandom(seed,p1,p2);
			}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	// ------------------------------------------------------------------
	// ESERCIZIO 5 ----------------------------------------------------

	unsigned int M = 50000;
	unsigned int N = 100;
	unsigned int L = M/N;

	double *ave = new double[N];

	double stepsize = 1.15;

	double* base_pos = new double[3];
	base_pos[0] = 1;
	base_pos[1] = 1;
	base_pos[2] = 1;

	double base_record_pos[M][3];

	unsigned int k;
	unsigned int rejects = 0;
	double sum, r, last_r;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;

			last_r = sqrt(pow(base_pos[0],2)  + pow(base_pos[1],2) + pow(base_pos[2],2));
			
			base_pos = metropolis(base_state, base_pos, stepsize, rnd);

			r = sqrt(pow(base_pos[0],2)  + pow(base_pos[1],2) + pow(base_pos[2],2));

			if(r == last_r){
				rejects++;
			}

			base_record_pos[k][0] = base_pos[0];
			base_record_pos[k][1] = base_pos[1];
			base_record_pos[k][2] = base_pos[2];

			sum += r;

		}

		ave[i] = sum/L;
	}		

	cout << "Base State Acceptance Rate: " << (int)((1 - (double) rejects/M)*100) << "%" <<  endl;

	fstream fout;
	fout.open("base_record_pos.csv", fstream::out);
	fout << "x, y, z" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout << base_record_pos[i][0] << ", "
		<< base_record_pos[i][1] << ", "
		<< base_record_pos[i][2] << endl;
	}
	fout.close();

	data_blocking(N, ave, "base_state.csv");


	// -----------------------------------------------------------

	double excited_record_pos[M][3];

	N = 200;
	L = M/N;

	ave = new double[N];

	stepsize = 2.7;

	double* ex_pos = new double[3];
	ex_pos[0] = 1;
	ex_pos[1] = 1;
	ex_pos[2] = 4;

	rejects = 0;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;

			last_r = sqrt(pow(ex_pos[0],2)  + pow(ex_pos[1],2) + pow(ex_pos[2],2));
			
			ex_pos = metropolis(excited_state, ex_pos, stepsize, rnd);

			r = sqrt(pow(ex_pos[0],2)  + pow(ex_pos[1],2) + pow(ex_pos[2],2));

			if(r == last_r){
				rejects++;
			}

			excited_record_pos[k][0] = ex_pos[0];
			excited_record_pos[k][1] = ex_pos[1];
			excited_record_pos[k][2] = ex_pos[2];

			sum += r;

		}

		ave[i] = sum/L;
	}		

	cout << "Excited State Acceptance Rate: " << (int)((1 - (double) rejects/M)*100) << "%" <<  endl;


	fout.open("excited_record_pos.csv", fstream::out);
	fout << "x, y, z" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout << excited_record_pos[i][0] << ", "
		<< excited_record_pos[i][1] << ", "
		<< excited_record_pos[i][2] << endl;
	}
	fout.close();

	data_blocking(N, ave, "excited_state.csv");

	return 0;
}

double error(double av, double av2, int n){
	if (n==0) return 0;
	else return sqrt((av2 - pow(av,2))/n);
};

void data_blocking(unsigned int N, double* ave, string fname){

	double * sum_prog = new double[N];
	double * su2_prog = new double[N];
	double * err_prog = new double[N];

	double av;
	for(unsigned int i = 0; i<N; i++){
		for(unsigned int j=0; j<i+1; j++){
			av = ave[j];
			sum_prog[i] += av;
			su2_prog[i] += pow(av,2);
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}

	fstream fout;
	fout.open(fname, fstream::out);
	fout << "block, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << i+1 << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	fout.close();

	delete sum_prog;
	delete su2_prog;
	delete err_prog;

}



double base_state(double x, double y, double z){

	double a0 = 1;

	double r = sqrt(x*x + y*y + z*z);

	double wave = pow(a0, -3/2)*exp(-a0*r)/sqrt(M_PI);

	return pow(wave,2);
}

/*
// D excited state
double excited_state(double x, double y, double z){

	double a0 = 1;

	double r = sqrt(x*x + y*y + z*z);
	double theta = acos(z/r);

	double wave = pow(a0, -3/2)*sqrt(2)/81/sqrt(M_PI)*(6-r/a0)*r/a0*exp(-r/3/a0)*cos(theta);

	return pow(wave,2);
}
*/

// P excited state
double excited_state(double x, double y, double z){

	double a0 = 1;

	double r = sqrt(x*x + y*y + z*z);
	double theta = acos(z/r);

	double wave = pow(a0, -5/2)*r*exp(-a0*r/2)*cos(theta)/8*sqrt(2/M_PI);

	return pow(wave,2);
}


double* metropolis(double (wave)(double, double, double), double* old_pos, double a, Random *rnd){

	double * new_pos = new double[3];

	double x = old_pos[0];
	double y = old_pos[1];
	double z = old_pos[2];

	double theta = rnd->Angle()/2;
	double phi = rnd->Angle();

	double nx = x + abs(a)*sin(theta)*cos(phi);
	double ny = y + abs(a)*sin(theta)*sin(phi);
	double nz = z + abs(a)*cos(theta);

	double A = min(1., wave(nx, ny, nz)/wave(x,y,z));

	if (rnd->Rannyu()<A){
		new_pos[0] = nx;
		new_pos[1] = ny;
		new_pos[2] = nz;
	} else {
		new_pos[0] = x;
		new_pos[1] = y;
		new_pos[2] = z;
	}

	return new_pos;

}


/****************************************************************
*****************************************************************
		_/		_/	_/_/_/	_/			 Numerical Simulation Laboratory
	 _/_/  _/ _/			 _/				Physics Department
	_/	_/_/		_/		_/			 Universita' degli Studi di Milano
 _/		 _/				_/ _/				Prof. D.E. Galli
_/		_/	_/_/_/	_/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/