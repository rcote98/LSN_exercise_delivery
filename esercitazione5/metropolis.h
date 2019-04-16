#ifndef __Metropolis__
#define __Metropolis__

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "random.h"

using namespace std;

class Metropolis3D{

public:

    Metropolis3D();
    Metropolis3D(double *start_pos);
    ~Metropolis3D(){};

    void Reset();
    double CalibrateFixed(double estimate ,double (func)(double, double, double));
    double CalibrateGaussian(double estimate ,double (func)(double, double, double));

    void GaussianStep(double sigma, double (func)(double, double, double));
    void FixedStep(double stepsize, double (func)(double, double, double));

    double* GetCurrentPos() const;
    double GetAcceptanceRate() const;

private:

    Random *rnd;
    int steps, rejects;
    double* initial_pos;
    double* current_pos;

};

#endif // __Metropolis__