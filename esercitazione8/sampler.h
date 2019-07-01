#ifndef __WSampler__
#define __WSampler__

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "random.h"
#include "metropolis.h"

using namespace std;

class WSampler{

public:

    WSampler();
    ~WSampler(){};

    double Measure();
    void Measure(string fname);
    void WaveHist(string fname, unsigned int n);

    static double Potential(double x);
    static double Wave(double x, double mu, double sigma);
    static double WaveSecond(double x, double mu, double sigma);
    
    static double Density(vector<double> pos, vector<double> params);
    static double Hamiltonian(vector<double> pos, vector<double> params);

    void SetMu(double nmu){mu = nmu;};
    void SetSigma(double nsigma){sigma = nsigma;};

    double GetMu() const {return mu;};
    double GetSigma() const {return sigma;};

    double error(double av, double av2, int n);

private:

    Metropolis* met;

    const unsigned int M = 50000;

    const double START = 1;
    vector<double> pos;

    double mu, sigma;
    double stepsize;
    double ham;

};

#endif // __WSampler__