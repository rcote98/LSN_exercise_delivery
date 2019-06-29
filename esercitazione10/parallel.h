#ifndef __TSProblem__
#define __TSProblem__

#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <string>

#include <cmath>
#include <vector>
#include <algorithm>

#include "random.h"

using namespace std;

class Annealing{

public:

    Annealing(vector< vector<double> > c, double temp_max, double temp_min, unsigned int iterations);
    ~Annealing(){};

    // temperature scheduler

    double LinearScheduler(unsigned int iteration);

    // iteration stuff
    void Move();

    // path generation
    vector<unsigned int> GeneratePath();

    // path mutation
    vector<unsigned int> SwapMutation(vector<unsigned int> path);

    // safety checks
    bool CheckPath();
    void PrintPath();
    bool CheckCity(unsigned int i);

    // aux functions
    double CostCalc(vector<unsigned int> path);
    double Norm(unsigned int i, unsigned int j);

    // getters
    vector<unsigned int> GetPath() const {return path;};
    unsigned int GetIteration() const {return iteration;};
    double GetCost() const {return curr_cost;};
    vector< vector<double> > GetCities() const {return cities;};

private:

    unsigned const int INIT_MUT = 60;

    unsigned int iters;
    unsigned int iteration;
    unsigned int ncity;

    Random *rnd;

    double tmax, tmin;
    double curr_cost;
    double beta;


    vector< unsigned int > path; 
    vector< vector< double > > cities;
    
};

#endif 