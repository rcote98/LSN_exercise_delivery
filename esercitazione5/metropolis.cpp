#include "metropolis.h"

Metropolis3D :: Metropolis3D(){

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

    current_pos = new double[3];
	initial_pos = new double[3];
	

    for(unsigned int i = 0; i < 3; i++){
        current_pos[i] = 0;
		initial_pos[i] = 0;
	}

}

Metropolis3D :: Metropolis3D(double *start_pos){

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

    current_pos = new double[3];
	initial_pos = new double[3];
	

    for(unsigned int i = 0; i < 3; i++){
        current_pos[i] = start_pos[i];
		initial_pos[i] = start_pos[i];
	}

}

void Metropolis3D :: Reset(){

	steps = 0;
	rejects = 0;

	for(unsigned int i = 0; i < 3; i++){
        current_pos[i] = initial_pos[i];
	}

}

double Metropolis3D :: CalibrateFixed(double estimate ,double (func)(double, double, double)){

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

double Metropolis3D :: CalibrateGaussian(double estimate ,double (func)(double, double, double)){

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

void Metropolis3D :: GaussianStep(double sigma, double (func)(double, double, double)){

	double x = current_pos[0];
	double y = current_pos[1];
	double z = current_pos[2];

	double nx = x + rnd->Gauss(0, sigma);
	double ny = y + rnd->Gauss(0, sigma);
	double nz = z + rnd->Gauss(0, sigma);

	double A = min(1., func(nx, ny, nz)/func(x,y,z));

	if (rnd->Rannyu()<A){
		current_pos[0] = nx;
		current_pos[1] = ny;
		current_pos[2] = nz;
	} else {
		rejects ++;
	}

	steps++;

}

void Metropolis3D :: FixedStep(double stepsize, double (func)(double, double, double)){

	double x = current_pos[0];
	double y = current_pos[1];
	double z = current_pos[2];

	double theta = rnd->Angle()/2;
	double phi = rnd->Angle();

	double nx = x + abs(stepsize)*sin(theta)*cos(phi);
	double ny = y + abs(stepsize)*sin(theta)*sin(phi);
	double nz = z + abs(stepsize)*cos(theta);

	double A = min(1., func(nx, ny, nz)/func(x,y,z));

	if (rnd->Rannyu()<A){
		current_pos[0] = nx;
		current_pos[1] = ny;
		current_pos[2] = nz;
	} else {
		rejects ++;
	}

	steps++;

}

double* Metropolis3D :: GetCurrentPos() const {

    double * return_pos = new double[3];

    for(unsigned int i = 0; i < 3; i++){
        return_pos[i] = current_pos[i];
    }

    return return_pos;

}

double Metropolis3D :: GetAcceptanceRate() const {

    if(steps == 0){
            return 0;
    }

    return (double)(steps - rejects)/(double)steps;

}


