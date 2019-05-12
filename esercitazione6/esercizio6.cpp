/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

#include "esercizio6.h"


int main(int argc, char *argv[])
{ 
    cli_argc = argc;
    cli_argv = argv;
    verbose = false;

    Input(); //Inizialization

    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep){
            Move(metro);
            Measure();
            Accumulate(); //Update block averages
        }

        Averages(iblk);   //Print results for current block
  
    }
    
    ConfFinal(); //Write final configuration

  return 0;
}

void Input(void) {

    ifstream ReadInput;

    if(verbose){
        cout << "Classic 1D Ising model             " << endl;
        cout << "Monte Carlo simulation             " << endl << endl;
        cout << "Nearest neighbour interaction      " << endl << endl;
        cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
        cout << "The program uses k_B=1 and mu_B=1 units " << endl;
    }
 

    // RNG SETUP ------------------------r-------------------------------

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
    
    //Read input informations
    
    ReadInput.open("input.dat");

    ReadInput >> temp;
    ReadInput >> nspin;
    ReadInput >> J;
    ReadInput >> h;
    ReadInput >> metro; // if=1 Metropolis else Gibbs
    ReadInput >> nblk;
    ReadInput >> nstep;
    ReadInput >> RESTART;

    if( cli_argc == 2){
        temp = atof(cli_argv[1]);
    }

    if( cli_argc == 3){
        temp = atof(cli_argv[1]);
        h = atof(cli_argv[2]);
    }

    bet = 1.0/temp;

    if(metro==1) cout << "The program performs Metropolis moves" << endl;
    else cout << "The program performs Gibbs moves" << endl;
    cout << "Temperature = " << temp << endl;
    cout << "Number of spins = " << nspin << endl;
    cout << "Exchange interaction = " << J << endl;
    cout << "External field = " << h << endl << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl;
    

    ReadInput.close();

    accepted = 0;
    attempted = 0;

    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility
 
    n_props = 4; //Number of observables

    //initial configuration
    if(RESTART){
        for (int i=0; i<nspin; ++i)
        {
                if(rnd->Rannyu() >= 0.5) s[i] = 1;
                else s[i] = -1;
        }
    } else
    {
        ReadInput.open("config.final");
        for (int i=0; i<nspin; ++i)
        {
            ReadInput >> s[i];
        }
        ReadInput.close();
    }
    
    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl << endl;
}


void Move(int metro){

    int o;
    double p, r, energy_old, energy_new;
    double energy_up, energy_down;

    for(int i=0; i<nspin; ++i)
    {
        //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        o = (int)(rnd->Rannyu()*nspin);

        if(metro==1) //Metropolis flipper
        {
        
            energy_old = Boltzmann(s[o], o);
            energy_new = Boltzmann(-1*s[o], o);

            if (energy_new - energy_old < 0){

                s[o] = -1*s[o];
                accepted++;
            
            } else {

                p = exp(-bet*(energy_new - energy_old));
                r = rnd->Rannyu();

                if (r < p) {
                    s[o] = -1*s[o];
                    accepted++;
                }
            }

            attempted++;

        }
        else //Gibbs sampling
        {
            
            energy_up = Boltzmann(1, o);
            energy_down = Boltzmann(-1, o);

            p = 1./(1.+ exp(-bet*(energy_down-energy_up)));
            r = rnd->Rannyu();


            if (r < p){

                s[o] = 1;
                accepted++;
            
            } else {
                
                s[o] = -1;

            }

            attempted++;
        }
    }
}


double Boltzmann(int sm, int ip)
{
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

void Measure()
{
    double u = 0.0; // energy
    double u2 = 0.0; // energy squared 
    double si = 0;

    //cycle over spins
    for (int i=0; i<nspin; ++i)
    {
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        u2 += pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]),2);
        si += s[i];
    }
    
    walker[iu] = u;
    walker[ic] = u2;
    walker[im] = si;
    walker[ix] = bet*si*si;
    
}


void Reset(int iblk) //Reset block averages
{
   
    if(iblk == 1)
    {
        for(int i=0; i<n_props; ++i)
        {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for(int i=0; i<n_props; ++i)
    {
        blk_av[i] = 0;
    }

    blk_norm = 0;
    attempted = 0;
    accepted = 0;

}



void Accumulate(void) //Update block averages
{
    for(int i=0; i<n_props; ++i)
    {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
        
    ofstream Ene, Heat, Mag, Chi;
    const int wd=15;
    

    cout << "Block number:  " << setw(2) << iblk << "; Acceptance rate: " << accepted/attempted << endl;
    

    //Energy
    stima_u = blk_av[iu]/blk_norm/(double)nspin; 

    glob_av[iu]  += stima_u;
    glob_av2[iu] += pow(stima_u, 2);
    err_u=error(glob_av[iu],glob_av2[iu],iblk);


    Ene.open("output.ene.0",ios::app);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    //Heat
    stima_c = bet*bet*(blk_av[ic]/blk_norm/(double)nspin - pow(blk_av[iu]/blk_norm/(double)nspin,2)); 

    glob_av[ic]  += stima_c;
    glob_av2[ic] += pow(stima_c, 2);
    err_c=error(glob_av[ic],glob_av2[ic],iblk);


    Heat.open("output.heat.0",ios::app);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    //Mag
    stima_m = blk_av[im]/blk_norm/(double)nspin; 

    glob_av[im]  += stima_m;
    glob_av2[im] += pow(stima_m, 2);
    err_m=error(glob_av[im],glob_av2[im],iblk);


    Mag.open("output.mag.0",ios::app);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    //Chi
    stima_x = blk_av[ix]/blk_norm/(double)nspin; 

    glob_av[ix]  += stima_x;
    glob_av2[ix] += pow(stima_x, 2);
    err_x=error(glob_av[ix],glob_av2[ix],iblk);


    Chi.open("output.chi.0",ios::app);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();


}

void ConfFinal(void)
{
    ofstream WriteConf;

    cout << endl << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    for (int i=0; i<nspin; ++i)
    {
        WriteConf << s[i] << endl;
    }
    WriteConf.close();

    rnd->SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
