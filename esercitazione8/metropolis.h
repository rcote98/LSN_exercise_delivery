#ifndef __Metropolis__
#define __Metropolis__

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "random.h"

using namespace std;

class Metropolis{

public:

    Metropolis(int dim);
    Metropolis(vector <double> start_pos);
    ~Metropolis(){};

    void Reset();
    double CalibrateFixed(double estimate ,double (func)(vector <double> pos, vector <double> params));
    double CalibrateGaussian(double estimate ,double (func)(vector <double> pos, vector <double> params));

    void GaussianStep(double sigma, double (func)(vector <double> pos, vector <double> params));
    void FixedStep(double stepsize, double (func)(vector <double> pos, vector <double> params));

    void SetParams(vector <double> p){params = p;};

    vector <double> GetCurrentPos() const;
    double GetAcceptanceRate() const;

private:

    Random *rnd;
    
    unsigned int dim;
    int steps, rejects;

    vector <double> initial_pos;
    vector <double> current_pos;

    vector <double> params; // parameters that dont vary with steps
};

#endif // __Metropolis__