/****************************************************************
*****************************************************************
		_/		_/	_/_/_/	_/			 Numerical Simulation Laboratory
	 _/_/  _/ _/			 _/				Physics Department
	_/	_/_/		_/		_/			 Universita' degli Studi di Milano
 _/		 _/				_/ _/				Prof. D.E. Galli
_/		_/	_/_/_/	_/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>			// cin, cout: Standard Input/Output Streams Library
#include <fstream>			// Stream class to both read and write from/to files.
#include <cmath>            // Math stuff

#include "random.h"
#include "metropolis.h"

using namespace std;

double error(double av, double av2, int n);
void data_blocking(unsigned int N, double* ave, string fname);

double base_state(double x, double y, double z);
double P_excited_state(double x, double y, double z);
double D_excited_state(double x, double y, double z);

int main(){
 
	// ##### BASE STATE #####

	// Blocking method stuff
	unsigned int M = 50000;
	unsigned int N = 100;
	unsigned int L = M/N;

	double *ave = new double[N];

	// Starting position
	double* pos = new double[3];
	pos[0] = 0;
	pos[1] = 0;
	pos[2] = 1.5;
	Metropolis3D base(pos); 

	double step_estimate = 1.1;
	double sigma_estimate = 1.1;

	// Metropolis Parameters
	double stepsize = base.CalibrateFixed(step_estimate, base_state);
	double sigma = base.CalibrateGaussian(sigma_estimate, base_state);
	double fixed_accept;
	double gaussian_accept;

	// Position Vectors
	double posisiton_record[M][3];

	// Fixed Step Calculation
	unsigned int k;
	double sum, r;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;
			
			base.FixedStep(stepsize, base_state);
			pos = base.GetCurrentPos();

			r = sqrt(pow(pos[0],2)  + pow(pos[1],2) + pow(pos[2],2));

			posisiton_record[k][0] = pos[0];
			posisiton_record[k][1] = pos[1];
			posisiton_record[k][2] = pos[2];

			sum += r;

		}
		ave[i] = sum/L;
	}	

	fixed_accept = base.GetAcceptanceRate();
	data_blocking(N, ave, "fixed_base_state.csv");
	base.Reset();

	// Gaussian Step Calculation
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;
			
			base.GaussianStep(sigma, base_state);
			pos = base.GetCurrentPos();

			r = sqrt(pow(pos[0],2)  + pow(pos[1],2) + pow(pos[2],2));

			sum += r;

		}
		ave[i] = sum/L;
	}	

	gaussian_accept = base.GetAcceptanceRate();
	data_blocking(N, ave, "gaussian_base_state.csv");

	// Console Output
	cout << "#### BASE STATE ####" << endl;
	cout << "Step     = " << stepsize << endl;
	cout << "Sigma    = " << sigma << endl;
	cout << "Fixed AR = " << (int) (fixed_accept*100) << "%" << endl;
	cout << "Gauss AR = " << (int) (gaussian_accept*100) << "%" << endl << endl;

	// Position File Output
	fstream fout;
	fout.open("base_record_pos.csv", fstream::out);
	fout << "x, y, z" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout << posisiton_record[i][0] << ", " << posisiton_record[i][1] << ", " << posisiton_record[i][2] << endl;
	}
	fout.close();



	// ----------- ##### P EXCITED STATE ##### -----------

	N = 200;
	L = M/N;

	ave = new double[N];

	// Starting position
	pos[0] = 0;
	pos[1] = 0;
	pos[2] = 5;
	Metropolis3D P_excited(pos); 

	step_estimate = 2.5;
	sigma_estimate = 2.1;


	// Metropolis Parameters
	stepsize = P_excited.CalibrateFixed(step_estimate, P_excited_state);
	sigma = P_excited.CalibrateGaussian(sigma_estimate, P_excited_state);

	// Fixed Step Calculation
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;
			
			P_excited.FixedStep(stepsize, P_excited_state);
			pos = P_excited.GetCurrentPos();

			r = sqrt(pow(pos[0],2)  + pow(pos[1],2) + pow(pos[2],2));

			posisiton_record[k][0] = pos[0];
			posisiton_record[k][1] = pos[1];
			posisiton_record[k][2] = pos[2];

			sum += r;

		}
		ave[i] = sum/L;
	}	

	fixed_accept = P_excited.GetAcceptanceRate();
	data_blocking(N, ave, "fixed_p_excited_state.csv");
	P_excited.Reset();

	// Gaussian Step Calculation
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;
			
			P_excited.GaussianStep(sigma, P_excited_state);
			pos = P_excited.GetCurrentPos();

			r = sqrt(pow(pos[0],2)  + pow(pos[1],2) + pow(pos[2],2));

			sum += r;

		}
		ave[i] = sum/L;
	}	

	gaussian_accept = P_excited.GetAcceptanceRate();
	data_blocking(N, ave, "gaussian_p_excited_state.csv");

	// Console Output
	cout << "#### P EXCITED STATE ####" << endl;
	cout << "Step     = " << stepsize << endl;
	cout << "Sigma    = " << sigma << endl;
	cout << "Fixed AR = " << (int) (fixed_accept*100) << "%" << endl;
	cout << "Gauss AR = " << (int) (gaussian_accept*100) << "%" << endl << endl;

	// Position File Output
	fout.open("p_excited_record_pos.csv", fstream::out);
	fout << "x, y, z" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout << posisiton_record[i][0] << ", " << posisiton_record[i][1] << ", " << posisiton_record[i][2] << endl;
	}
	fout.close();


	// ----------- ##### D EXCITED STATE ##### -----------

	N = 200;
	L = M/N;

	ave = new double[N];

	// Starting position
	pos[0] = 0;
	pos[1] = 0;
	pos[2] = 12.5;
	Metropolis3D D_excited(pos); 

	step_estimate = 3.5;
	sigma_estimate = 3.5;

	// Metropolis Parameters
	stepsize = D_excited.CalibrateFixed(step_estimate, D_excited_state);
	sigma = D_excited.CalibrateGaussian(sigma_estimate, D_excited_state);

	// Fixed Step Calculation
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;
			
			D_excited.FixedStep(stepsize, D_excited_state);
			pos = D_excited.GetCurrentPos();

			r = sqrt(pow(pos[0],2)  + pow(pos[1],2) + pow(pos[2],2));

			posisiton_record[k][0] = pos[0];
			posisiton_record[k][1] = pos[1];
			posisiton_record[k][2] = pos[2];

			sum += r;

		}
		ave[i] = sum/L;
	}	

	fixed_accept = D_excited.GetAcceptanceRate();
	data_blocking(N, ave, "fixed_d_excited_state.csv");
	D_excited.Reset();

	// Gaussian Step Calculation
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			k = i*L + j;
			
			D_excited.GaussianStep(sigma, D_excited_state);
			pos = D_excited.GetCurrentPos();

			r = sqrt(pow(pos[0],2)  + pow(pos[1],2) + pow(pos[2],2));

			sum += r;

		}
		ave[i] = sum/L;
	}	

	gaussian_accept = D_excited.GetAcceptanceRate();
	data_blocking(N, ave, "gaussian_d_excited_state.csv");

	// Console Output
	cout << "#### D EXCITED STATE ####" << endl;
	cout << "Step     = " << stepsize << endl;
	cout << "Sigma    = " << sigma << endl;
	cout << "Fixed AR = " << (int) (fixed_accept*100) << "%" << endl;
	cout << "Gauss AR = " << (int) (gaussian_accept*100) << "%" << endl << endl;

	// Position File Output
	fout.open("d_excited_record_pos.csv", fstream::out);
	fout << "x, y, z" << endl;
	for(unsigned int i = 0; i < M; i++){
		fout << posisiton_record[i][0] << ", " << posisiton_record[i][1] << ", " << posisiton_record[i][2] << endl;
	}
	fout.close();



	// ------------------------------------------------------------------------
	
	return 0;

}

// ###########################################################################3

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

// P excited state
double P_excited_state(double x, double y, double z){

	double a0 = 1;

	double r = sqrt(x*x + y*y + z*z);
	double theta = acos(z/r);

	double wave = pow(a0, -5/2)*r*exp(-a0*r/2)*cos(theta)/8*sqrt(2/M_PI);

	return pow(wave,2);
}

// D excited state
double D_excited_state(double x, double y, double z){

	double a0 = 1;

	double r = sqrt(x*x + y*y + z*z);
	double theta = acos(z/r);

	double wave = pow(a0, -3/2)*sqrt(2)/81/sqrt(M_PI)*(6-r/a0)*r/a0*exp(-r/3/a0)*cos(theta);

	return pow(wave,2);
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