#include "metropolis.h"

Metropolis :: Metropolis(int dims){

    // RNG SETUP -------------------------------------------------------

	rnd = new Random();

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

    steps = 0;
    rejects = 0;
	dim = dims;

    current_pos.reserve(dim);
	initial_pos.reserve(dim);

	for(unsigned int i = 0; i < dim; i++){
		current_pos[i] = 0;
		initial_pos[i] = 0;
	}

}

Metropolis :: Metropolis(vector <double> start_pos){

    // RNG SETUP -------------------------------------------------------

	rnd = new Random();

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

    steps = 0;
    rejects = 0;
	dim = start_pos.size();

	current_pos = start_pos;
    initial_pos = start_pos;

}

void Metropolis :: Reset(){

	steps = 0;
	rejects = 0;

    current_pos = initial_pos;
}

double Metropolis :: CalibrateFixed(double estimate ,double (func)(vector <double> pos, vector <double> params)){

	unsigned int trials = 1000;
	double var_rate = 0.05;
	Reset();
	
	for(unsigned int i = 0; i < trials; i++)
	{
		FixedStep(estimate, func);
	}
	
	while(GetAcceptanceRate() < .47 || GetAcceptanceRate() > 0.53){

		if(GetAcceptanceRate() > .5){
			estimate += estimate*var_rate;
		} else estimate -= estimate*var_rate;

		Reset();

		for(unsigned int i = 0; i < trials; i++){
			FixedStep(estimate, func);
		}

	};


	return estimate;

}

double Metropolis :: CalibrateGaussian(double estimate ,double (func)(vector <double> pos, vector <double> params)){

	unsigned int trials = 1000;
	double var_rate = 0.05;
	Reset();
	
	for(unsigned int i = 0; i < trials; i++)
	{
		GaussianStep(estimate, func);
	}
	
	while(GetAcceptanceRate() < .47 || GetAcceptanceRate() > 0.53){

		if(GetAcceptanceRate() > .5){
			estimate += estimate*var_rate;
		} else estimate -= estimate*var_rate;

		Reset();

		for(unsigned int i = 0; i < trials; i++){
			GaussianStep(estimate, func);
		}

	};

	Reset();

	return estimate;
	
}

void Metropolis :: GaussianStep(double sigma, double (func)(vector <double> pos, vector <double> params)){

	vector<double> nx = current_pos;
	for(unsigned int i = 0; i < dim; i++){
		nx[i] += rnd->Gauss(0, sigma);
	}

	double A = min(1., func(nx, params)/func(current_pos, params));

	if (rnd->Rannyu()<A){
		current_pos = nx;
	} else {
		rejects ++;
	}

	steps++;

}

void Metropolis :: FixedStep(double stepsize, double (func)(vector <double> pos, vector <double> params)){

	vector<double> nx = current_pos;

	for(unsigned int i = 0; i < dim; i++){
		
		nx[i] += stepsize*rnd->Rannyu() - stepsize/2;

	}

	double A = min(1., func(nx, params)/func(current_pos, params));

	if (rnd->Rannyu()<A){
		current_pos = nx;
	} else {
		rejects ++;
	}

	steps++;

}

vector<double> Metropolis :: GetCurrentPos() const {

    vector<double> return_pos = current_pos;

    return return_pos;

}

double Metropolis :: GetAcceptanceRate() const {

    if(steps == 0){
            return 0;
    }

    return (double)(steps - rejects)/(double)steps;

}


