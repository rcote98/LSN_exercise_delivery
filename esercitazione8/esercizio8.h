/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

unsigned int arg_count;
char ** args;

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables


double x;
double p;
double mu, sigma;

// simulation

const int m_props=1000;

double p_trans;
double delta = 3;
double nblk = 30;
double nstep = 10000;

const int wd = 15;

int mcstep;
bool verbose;

int n_props, ih;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_H, err_H;

unsigned int accepted, attempted;

// constants
const double pi = 3.1415927;
const double m = 1;
const double hbar = 1;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
double Error(double,double,int);


// physics functions
double Hamiltonian(double x, double mu, double sigma);

double Potential(double x);
double Wave(double x, double mu, double sigma);
double Wave_second(double x, double mu, double sigma);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
