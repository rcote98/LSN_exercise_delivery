#include <iostream>			// cin, cout: Standard Input/Output Streams Library
#include <fstream>			// Stream class to both read and write from/to files.
#include <iomanip>
#include <cmath>            // Math stuff

#include "sampler.h"

using namespace std;

int main(){

	const double MU_MIN = 0.5;
	const double MU_MAX = 1;

	const double SIGMA_MIN = 0.5;
	const double SIGMA_MAX = 1;

	double mu_range = MU_MAX - MU_MIN;
	double sg_range = SIGMA_MAX - SIGMA_MIN;

	unsigned int nmu = 30;
	unsigned int nsigma = 30;

	WSampler sampler = WSampler();

	double ham, mu, sigma;

	double bestham = 50;
	double bestmu = 0; 
	double bestsigma = 0;

	cout << "GROUND STATE OPTIMIZATION W/VARIATIONAL MONTE CARLO" << endl << endl;
	cout << "Total number of states to sample: " << nmu*nsigma << endl << endl;
	cout << "Beginning optimization process..." << endl;

	fstream fout;
	fout.open("optimizing.dat", fstream::out);

	fout << setw(20) << "mu" << setw(20) << "sigma" << setw(20) << "energy" << endl;
	for(unsigned int i = 0; i < nmu; i++){

		if(i== nmu/4) cout << "25% done..." << endl;
		if(i== nmu/2) cout << "50% done..." << endl;
		if(i== 3*nmu/4) cout << "75% done..." << endl;

		for (unsigned int j = 0; j < nsigma; j++){

			mu = MU_MIN + (double)i*mu_range/nmu;
			sigma = SIGMA_MIN + (double)j*sg_range/nsigma;

			sampler.SetMu(mu);
			sampler.SetSigma(sigma);

			ham = sampler.Measure();

			//cout << setw(20) << mu << setw(20) << sigma << setw(20) << ham << endl;
			fout << setw(20) << mu << setw(20) << sigma << setw(20) << ham << endl;

			if(ham < bestham){
				bestham = ham;
				bestmu = mu;
				bestsigma = sigma;
			}

		}
	}

	fout.close();

	cout << "DONE!" << endl << endl;

	cout << "Best Mu: " << bestmu << endl;
	cout << "Best Sigma: " << bestsigma << endl;
	cout << "Energy: " << bestham << endl << endl;

	cout << "Sampling WF with best parameters..." << endl;
	sampler.SetMu(bestmu);
	sampler.SetSigma(bestsigma);
	cout << "DONE!" << endl;


	sampler.Measure("hamiltonian.dat");
	sampler.WaveHist("wave_hist.dat", 200000);
	
	return 0;

}