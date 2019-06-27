/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <stdlib.h>
#include <iostream>  
#include <fstream>      
#include <cmath>
#include <string>

#include "random.h"

using namespace std;

unsigned int arg_count;
char ** args;

// random numbers

Random * rnd;

// parameters, observables
//const int m_props=4;
double stima_pot, stima_kin, stima_etot;
double stima_temp, stima_pres, stima_gofr;

// averages
double acc,att;

// g(r)
unsigned const int nbins = 100;
double gofr_hist[nbins];
double interval_min, interval_max; // bin max/min
double norm; // histogram normalization
double max_radius;


//configuration
unsigned const int m_part=108;
double x[m_part],y[m_part],z[m_part];
double xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
unsigned int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
unsigned int nstep, nblocks, iprint;
unsigned int measure_step=10;
double delta;
bool verbose;
bool RESTART;

// functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Averages();
void Measure(void);

double Force(unsigned int, unsigned int);
double Pbc(double);

void data_blocking(unsigned int N, double* ave, string fname);
double error(double av, double av2, int n);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
