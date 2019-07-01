#include "sampler.h"


WSampler :: WSampler(){

	mu = 0;
	sigma = 0;
	ham = 0;

	vector<double> start_pos = {START};

	met = new Metropolis(start_pos);
}


double WSampler::Measure(){
	
	met->Reset();

	// Blocking method stuff
	unsigned int N = 100;
	unsigned int L = M/N;

	double *ave = new double[N];

	vector<double> pars = {mu, sigma};
	met->SetParams(pars);

	stepsize = met->CalibrateFixed(mu, Density);

	double h, sum;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){

			met->FixedStep(stepsize, Density);
			pos = met->GetCurrentPos();

			h = Hamiltonian(pos, pars);

			sum += h;

		}
		ave[i] = sum/L;
	}	

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


	double best_stime = sum_prog[N-1];

	delete sum_prog;
	delete su2_prog;
	delete err_prog;

	met->Reset();

	return best_stime;

}


void WSampler::Measure(string fname){

	met->Reset();

	// Blocking method stuff
	unsigned int M = 50000;
	unsigned int N = 100;
	unsigned int L = M/N;

	double *ave = new double[N];

	vector<double> pars = {mu, sigma};
	met->SetParams(pars);

	stepsize = met->CalibrateFixed(mu, Density);

	double h, sum;
	for(unsigned int i = 0; i<N; i++){
		sum = 0;

		for(unsigned int j = 0; j<L; j++){
			
			met->FixedStep(stepsize, Density);
			pos = met->GetCurrentPos();

			h = Hamiltonian(pos, pars);

			sum += h;

		}
		ave[i] = sum/L;
	}	

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
	fout << "block, ave, err" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << i+1 << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	fout.close();

	delete sum_prog;
	delete su2_prog;
	delete err_prog;

	met->Reset();


}

void WSampler::WaveHist(string fname, unsigned int n){

	met->Reset();

	vector<double> pars = {mu, sigma};
	met->SetParams(pars);

	stepsize = met->CalibrateFixed(mu, Density);

	fstream fout;
	fout.open(fname, fstream::out);

	for(unsigned int i = 0; i<n; i++){

		met->FixedStep(stepsize, Density);
		pos = met->GetCurrentPos();

		fout << pos[0] << endl;
	}	
	
	fout.close();

	met->Reset();


}

double WSampler::Potential(double x){

	//return 0;

	return pow(x,4) - pow(x,2)*5/2;

}


double WSampler::Wave(double x, double mu, double sigma){

	return exp(-pow(x-mu,2)/(2*sigma*sigma)) + exp(-pow(x+mu,2)/(2*sigma*sigma));

}


double WSampler::WaveSecond(double x, double mu, double sigma){

	return -Wave(x, mu, sigma)/(sigma*sigma)
	+ exp(-pow(x+mu,2)/(2*sigma*sigma))*pow((x+mu)/(sigma*sigma),2)
	+ exp(-pow(x-mu,2)/(2*sigma*sigma))*pow((x-mu)/(sigma*sigma),2);

}


double WSampler::Density(vector<double> pos, vector<double> params){

	double x = pos[0];
	double mu = params[0];
	double sigma = params[1];

	return Wave(x, mu, sigma)*Wave(x, mu, sigma);

}


double WSampler::Hamiltonian(vector<double> pos, vector<double> params){

	double hbar = 1;
	double m = 1;

	double x = pos[0];
	double mu = params[0];
	double sigma = params[1];

	return (-0.5*hbar*hbar*WaveSecond(x, mu, sigma)/m + Potential(x)*Wave(x, mu, sigma))/Wave(x, mu, sigma); 

}


double WSampler::error(double av, double av2, int n){
	if (n==0) return 0;
	else return sqrt((av2 - pow(av,2))/n);
}