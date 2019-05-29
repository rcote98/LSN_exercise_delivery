/****************************************************************
*****************************************************************
		_/		_/	_/_/_/	_/			 Numerical Simulation Laboratory
	 _/_/  _/ _/			 _/				Physics Department
	_/	_/_/		_/		_/			 Universita' degli Studi di Milano
 _/		 _/				_/ _/				Prof. D.E. Galli
_/		_/	_/_/_/	_/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>

#include "esercizio8.h"
#include "random.h"


using namespace std;

int main(int argc, char * argv[]){

	verbose = false;

	arg_count = argc;
	args = argv;

	Input();

	for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
	{
		Reset(iblk); //Reset block averages
		for(int istep=1; istep <= nstep; ++istep)
		{

			Move();
			Measure();
			Accumulate(); //Update block averages

		}
		Averages(iblk);   //Print results for current block
	}
}

void Input(void)
{



	if(verbose){

		cout << "lul" << endl;

	}
	
	//Read seed for random numbers
	 int p1, p2;
	 ifstream Primes("Primes");
	 Primes >> p1 >> p2 ;
	 Primes.close();

	 ifstream input("seed.in");
	 input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	 rnd.SetRandom(seed,p1,p2);
	 input.close();
	
	
	// starting positions

	x = 1;

	mu = 6;
	sigma = 2;


	if (arg_count == 3){

		mu = atof(args[1]);
		sigma = atof(args[2]);

	}

	cout << "mu: " << mu << endl;
	cout << "sigma: " << sigma << endl;

	p_trans = 0.5;

	p = Wave(x, mu, sigma)*Wave(x, mu, sigma);

	cout << "The program perform Metropolis moves with uniform translations" << endl;
	cout << "Moves parameter = " << delta << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	

	//Prepare arrays for measurements
	ih = 0; //Potential energy
	n_props = 1; //Number of observables

}


void Reset(int iblk) //Reset block averages
{

	 if(iblk == 1)
	 {
			 for(int i=0; i<n_props; ++i)
			 {
					 glob_av[i] = 0;
					 glob_av2[i] = 0;
			 }
	 }

	 for(int i=0; i<n_props; ++i)
	 {
		 blk_av[i] = 0;
	 }

	 blk_norm = 0;
	 attempted = 0;
	 accepted = 0;
}

void Accumulate(void) //Update block averages
{

	 for(int i=0; i<n_props; ++i)
	 {
		 blk_av[i] = blk_av[i] + walker[i];
	 }
	 blk_norm = blk_norm + 1.0;

}

void Move(void)
{

	double energy;
	double pnew;
	double xnew;
	double r;

	xnew = x + delta*(rnd.Rannyu() - 0.5);
	pnew = Wave(xnew, mu, sigma)*Wave(xnew, mu, sigma);

	if(pnew > p){
		x = xnew;
		p = pnew;
		accepted++;

	} else {

		r = rnd.Rannyu();

		if(r < p_trans){
			x = xnew;
			p = pnew;
			accepted++;
		}
	}

	attempted ++;

}


void Measure()
{

	stima_H = p*Hamiltonian(x, mu, sigma)/Wave(x, mu, sigma);

	walker[ih] = stima_H;

	ofstream Prob;
	Prob.open("output.inst_pos.0",ios::app);

	Prob << setw(wd) << mcstep <<  setw(wd) << x << endl;

	Prob.close();
}


void Averages(int iblk) //Print results for current block
{
	ofstream Ham;
		
	cout << "Block number " << iblk << endl;
	cout << "Acceptance rate " << (double)accepted/(double)attempted << endl << endl;
		
	Ham.open("output.ham.0",ios::app);

	stima_H = blk_av[ih]/blk_norm;
	glob_av[ih] += stima_H;
	glob_av2[ih] += stima_H*stima_H;
	err_H=Error(glob_av[ih],glob_av2[ih],iblk);

	Ham << setw(wd) << iblk <<  setw(wd) << stima_H << setw(wd) << glob_av[ih]/(double)iblk << setw(wd) << err_H << endl;

	cout << "----------------------------" << endl << endl;

	Ham.close();
}

double Error(double sum, double sum2, int iblk)
{
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


//###################################################################################


double Hamiltonian(double x, double mu, double sigma){

	return -0.5*hbar*hbar/m*Wave_second(x, mu, sigma) + Potential(x);

}

double Potential(double x){

	return pow(x, 4) - 5/2*pow(x,2);

}

double Wave(double x, double mu, double sigma){

	return exp(-pow( x - mu ,2)/(2*pow( sigma ,2))) +  exp(-pow( x + mu ,2)/(2*pow( sigma ,2)));

}

double Wave_second(double x, double mu, double sigma){

	double positive = pow( (mu + x)/(sigma*sigma) ,2)*exp(-pow( x + mu ,2)/(2*pow( sigma ,2)));
	double negative = pow( (mu - x)/(sigma*sigma) ,2)*exp(-pow( x - mu ,2)/(2*pow( sigma ,2)));
	double common = 1/pow(sigma, 2)*(exp(-pow( x + mu ,2)/(2*pow( sigma ,2))) - exp(-pow( x - mu ,2)/(2*pow( sigma ,2))));

	return positive + negative + common;
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