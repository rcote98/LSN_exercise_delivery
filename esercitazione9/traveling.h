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

class TSProblem{

public:

    TSProblem(vector< vector<double> > c);
    ~TSProblem(){};

    // generation stuff

    void AdvanceGeneration();

    // path generation

    vector<unsigned int> GeneratePath();

    // safety checks

    bool CheckCity(unsigned int i);
    bool CheckPath(unsigned int i);


    // fitness stuff

    void CostCalc();
    void SortPopulation();

    // mutation stuff

    void Crossover(unsigned int i, unsigned int j);
    vector<unsigned int> SwapMutation(vector<unsigned int> path);
    vector<unsigned int> SwapMutation(vector<unsigned int> pathm, unsigned int m);
    vector<unsigned int> ShiftMutation(vector<unsigned int> pathm, unsigned int m);

    // statistics stuff

    double AverageFitness();
    double HalfAverageFitness();
    double AveragePathLength();
    double HalfAveragePathLength();

    // aux functions

    int Partition(int low, int high); 
    void QuickSort(int low, int high);
    double Norm(unsigned int i, unsigned int j);

    // printers

    void ShowPops(unsigned int n);

    // getters

    vector<unsigned int> GetBestPath() const {return paths[0];};
    unsigned int GetGeneration() const {return generation;};
    vector< vector<double> > GetCities() const {return cities;};
    vector< vector<unsigned int> > GetPaths() const {return paths;};


private:

    unsigned const int INIT_MUT = 60;
    unsigned const int MUT_GEN = 4000;
    unsigned const int CROSSOVERS = 100;
    const double M_EXP = 0.2;    
    const double C_EXP = 2;

    const double PSWAP = 1;
    const double PSHIFT = 1;
    const double PCROSS = 0.5;    

    unsigned int generation;
    unsigned int pop_size;
    unsigned int ncity;

    Random *rnd;
    vector< double > costs;
    vector< double > fitness;
    vector< vector< double > > cities;
    vector< vector< unsigned int > > paths; 
};

#endif 