/****************************************************************
*****************************************************************
		_/		_/	_/_/_/	_/			 Numerical Simulation Laboratory
	 _/_/  _/ _/			 _/				Physics Department
	_/	_/_/		_/		_/			 Universita' degli Studi di Milano
 _/		 _/				_/ _/				Prof. D.E. Galli
_/		_/	_/_/_/	_/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>			// srand, rand: to generate random number
#include <iostream>			// cin, cout: Standard Input/Output Streams Library
#include <fstream>			// Stream class to both read and write from/to files.
#include <cmath>			// rint, pow
#include "random.h"
#include "MolDyn_NVE.h"

using namespace std;

int main(){
 
	Input();						 //Inizialization
	int nconf = 1;
	for(int istep=1; istep <= nstep; ++istep){
		 Move();					 //Move particles with Verlet algorithm
		 if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
		 if(istep%10 == 0){
				Measure();		 //Properties measurement
				ConfXYZ(nconf);  //Write actual configuration in XYZ format
				nconf += 1;
		 }
		 
	}
	StatCalc();                  //Stats calculation
	ConfFinal();				 //Write final configuration to restart

	return 0;
}


void Input(void){ //Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf;
	double ep, ek, pr, et, vir;

	cout << "Classic Lennard-Jones fluid				" << endl;
	cout << "Molecular dynamics simulation in NVE ensemble	" << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	seed = 1;		 //Set seed for random numbers
	srand(seed); //Initialize random number generator
	
	ReadInput.open("input.dat"); //Read input

	ReadInput >> temp;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> iprint;
	ReadInput >> nblocks;

	// Variables for the blocking method
	etot_ave = new double[nblocks];
	epot_ave = new double[nblocks];
	ekin_ave = new double[nblocks];
	temp_ave = new double[nblocks];
	pres_ave = new double[nblocks];

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
	
	ReadInput >> RESTART;
	ReadInput >> REWRITE;
	ReadInput >> sim_name;

	ReadInput.close();

	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	ip = 4; //Pressure
	n_props = 5; //Number of observables

	if(RESTART){

		sim_name = "0"; // For file naming purposes

		// ----------------------------------------------------------------
		// RNG SETUP
		// ----------------------------------------------------------------

		Random rnd;

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
					rnd.SetRandom(seed,p1,p2);
				}
			}
		input.close();
		} else cerr << "PROBLEM: Unable to open seed.in" << endl;
		// ------------------------------------------------------------------


		//Read initial configuration
		cout << "Read initial configuration from file config/config.0 " << endl << endl;
		ReadConf.open("config/config.0");
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = x[i] * box;
			y[i] = y[i] * box;
			z[i] = z[i] * box;
		}
		ReadConf.close();


		//Prepare initial velocities
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i=0; i<npart; ++i){
			vx[i] = rnd.Rannyu() - 0.5;
			vy[i] = rnd.Rannyu() - 0.5;
			vz[i] = rnd.Rannyu() - 0.5;

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		double sumv2 = 0.0, fs;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);		// fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = x[i] - vx[i] * delta;
			yold[i] = y[i] - vy[i] * delta;
			zold[i] = z[i] - vz[i] * delta;
		}

	} else {

		string fname_actual = "config." + sim_name; 
		

		//Read initial position configuration
		cout << "Read initial 'position' from file" << fname_actual << endl << endl;
		ReadConf.open(fname_actual);
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = x[i] * box;
			y[i] = y[i] * box;
			z[i] = z[i] * box;
		}
		ReadConf.close();

		string fname_prev = "config.prev." + sim_name;

		//Read initial previous position configuration
		cout << "Read initial 'previous position' from file" << fname_prev << endl << endl;
		ReadConf.open(fname_prev);
		for (int i=0; i<npart; ++i){
			ReadConf >> xold[i] >> yold[i] >> zold[i];
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}
		ReadConf.close();


		//Prepare initial velocities
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		
		double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

		for(int j=0; j < 4; j++){ // repeat this process until temp is right

			for(int i=0; i<npart; ++i){ //Force acting on particle i
				fx[i] = Force(i,0);
				fy[i] = Force(i,1);
				fz[i] = Force(i,2);
			}

			for(int i=0; i<npart; ++i){ //Verlet integration scheme

				xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
				ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
				znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

				vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
				vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
				vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
			}



		}
		
		double sumv[3] = {0.0, 0.0, 0.0};


		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

		double sumv2 = 0.0, fs;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);		// fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = x[i] - vx[i] * delta;
			yold[i] = y[i] - vy[i] * delta;
			zold[i] = z[i] - vz[i] * delta;
		}


	}

	 return;
}


void Move(void){ //Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){ //Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){ //Verlet integration scheme

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}
	return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
	double f=0.0;
	double dvec[3], dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );	// distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );

			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
			}
		}
	}
	
	return f;
}

void Measure(){ //Properties measurement
	int bin;
	double v, t, vij;
	double dx, dy, dz, dr;
	ofstream Epot, Ekin, Etot, Temp;

	Epot.open("output_epot.dat",ios::app);
	Ekin.open("output_ekin.dat",ios::app);
	Temp.open("output_temp.dat",ios::app);
	Etot.open("output_etot.dat",ios::app);

	v = 0.0; //reset observables
	t = 0.0; 

	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){

		for (int j=i+1; j<npart; ++j){

		 dx = Pbc( x[i] - x[j] );
		 dy = Pbc( y[i] - y[j] );
		 dz = Pbc( z[i] - z[j] );

		 dr = dx*dx + dy*dy + dz*dz;
		 dr = sqrt(dr);

		 if(dr < rcut){
			 vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

			//Potential energy
			 v += vij;
		 }
		}					 
	}



	//Kinetic energy
	for (int i=0; i<npart; ++i){
	
		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	 
		stima_pot = v/(double)npart; //Potential energy
		stima_kin = t/(double)npart; //Kinetic energy
		stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
		stima_etot = (t+v)/(double)npart; //Total enery

		Epot << stima_pot  << endl;
		Ekin << stima_kin  << endl;
		Temp << stima_temp << endl;
		Etot << stima_etot << endl;

		Epot.close();
		Ekin.close();
		Temp.close();
		Etot.close();

		return;
	}
}

void StatCalc(void){


	ifstream Epot, Ekin, Etot, Temp;//, Pres;

	Epot.open("output_epot.dat");
	Ekin.open("output_ekin.dat");
	Temp.open("output_temp.dat");
	Etot.open("output_etot.dat");
	//Pres.open("output_pres.dat");

	unsigned int nmeasures = nstep/iprint;

	unsigned int k;
	double L = nmeasures/nblocks;
	double* sum = new double[4];


	double etot_meas[nmeasures];
	double ekin_meas[nmeasures];
	double epot_meas[nmeasures];
	double temp_meas[nmeasures];
	//double pres_meas[nmeasures];

	for(unsigned int i = 0; i < nmeasures; i++)
	{
			Epot >> epot_meas[i];
			Ekin >> ekin_meas[i];
			Temp >> temp_meas[i];
			Etot >> etot_meas[i];

			cout << epot_meas[i] << endl;
	}
	


	etot_ave = new double[nblocks];
	epot_ave = new double[nblocks];
	ekin_ave = new double[nblocks];
	temp_ave = new double[nblocks];
	pres_ave = new double[nblocks];

	for(unsigned int i = 0; i<nblocks; i++){
		
		for(unsigned int j = 0; j < 4; j++) sum[j] = 0;


		for(unsigned int j = 0; j < L; j++){

			k = j+1;

			sum[0] += epot_meas[k];
			sum[1] += ekin_meas[k];
			sum[2] += temp_meas[k];
			sum[3] += etot_meas[k];
		}

		epot_ave[i] = sum[0]/L;
		ekin_ave[i] = sum[1]/L;
		temp_ave[i] = sum[2]/L;
		etot_ave[i] = sum[3]/L;

	}

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();
	//Pres.close()

	
	unsigned int x[nblocks];
	for(unsigned int i = 0; i < nblocks; i++){
		x[i] = (i+1)*L;
	};
	

	ofstream Stats;

	double* ave;
	double* sum_prog = new double[nblocks];
	double* su2_prog = new double[nblocks];
	double* err_prog = new double[nblocks];

	// Total Energy #####################################################################
	ave = etot_ave;
	for(unsigned int i = 0; i<nblocks; i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += pow(ave[j],2);
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}

	Stats.open("ave_etot.out");
	Stats << "step, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < nblocks; i++){
		Stats << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	Stats.close();

	for(unsigned int i = 0; i < nblocks; i++){
		sum_prog[i]=0;
		su2_prog[i]=0;
		err_prog[i]=0;
	};

	// Potential Energy #####################################################################
	ave = epot_ave;
	for(unsigned int i = 0; i<nblocks; i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += pow(ave[j],2);
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}

	Stats.open("ave_epot.out");
	Stats << "step, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < nblocks; i++){
		Stats << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	Stats.close();

	for(unsigned int i = 0; i < nblocks; i++){
		sum_prog[i]=0;
		su2_prog[i]=0;
		err_prog[i]=0;
	};

	// Kinetic Energy #####################################################################
	ave = ekin_ave;
	for(unsigned int i = 0; i<nblocks; i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += pow(ave[j],2);
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}

	Stats.open("ave_ekin.out");
	Stats << "step, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < nblocks; i++){
		Stats << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	Stats.close();

	for(unsigned int i = 0; i < nblocks; i++){
		sum_prog[i]=0;
		su2_prog[i]=0;
		err_prog[i]=0;
	};

	// Temperature #####################################################################
	ave = temp_ave;
	for(unsigned int i = 0; i<nblocks; i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += pow(ave[j],2);
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}

	Stats.open("ave_temp.out");
	Stats << "step, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < nblocks; i++){
		Stats << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	Stats.close();

	for(unsigned int i = 0; i < nblocks; i++){
		sum_prog[i]=0;
		su2_prog[i]=0;
		err_prog[i]=0;
	};

	// Pressure #####################################################################
	ave = pres_ave;
	for(unsigned int i = 0; i<nblocks; i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += pow(ave[j],2);
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}

	Stats.open("ave_pres.out");
	Stats << "step, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < nblocks; i++){
		Stats << x[i] << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	Stats.close();

}


void ConfFinal(void){ //Write final configuration

	cout << endl;

	ofstream WriteConf;

	if(RESTART){

		cout << "Print final 'position' file config/config.final.0" << endl << endl;
		WriteConf.open("config/config.final.0");

		for (int i=0; i<npart; ++i){
			WriteConf << x[i]/box << "	 " <<  y[i]/box << "	 " << z[i]/box << endl;
		}
		WriteConf.close();

		cout << "Print final 'previous position' file config/config.final.prev.0" << endl << endl;
		WriteConf.open("config/config.final.prev.0");
		for (int i=0; i<npart; ++i){
			WriteConf << xold[i]/box << "	 " <<  yold[i]/box << "	 " << zold[i]/box << endl;
		}
		WriteConf.close();

	}else{

		if(REWRITE){
			string fname_final = "config/config." + sim_name;
			cout << "Overwriting initial 'position' file " << fname_final << endl << endl;
			
			WriteConf.open(fname_final);
			for (int i=0; i<npart; ++i){
				WriteConf << x[i]/box << "	 " <<  y[i]/box << "	 " << z[i]/box << endl;
			}
			WriteConf.close();

			string fname_final_prev = "config/config.prev." + sim_name;
			cout << "Overwriting initial 'previous position' file " << fname_final_prev << endl << endl;

			WriteConf.open(fname_final_prev);
			for (int i=0; i<npart; ++i){
				WriteConf << xold[i]/box << "	 " <<  yold[i]/box << "	 " << zold[i]/box << endl;
			}
			WriteConf.close();


		} else
		{
			string fname_final = "config.final." + sim_name;
			cout << "Print final 'position' file " << fname_final  << endl << endl;
			WriteConf.open(fname_final);

			for (int i=0; i<npart; ++i){
				WriteConf << x[i]/box << "	 " <<  y[i]/box << "	 " << z[i]/box << endl;
			}
			WriteConf.close();

			string fname_final_prev = "config/config.final.prev." + sim_name;
			cout << "Print final 'previous position' file " << fname_final_prev << endl << endl;

			WriteConf.open(fname_final_prev);
			for (int i=0; i<npart; ++i){
				WriteConf << xold[i]/box << "	 " <<  yold[i]/box << "	 " << zold[i]/box << endl;
			}
			WriteConf.close();
		}
		

	}
	


	WriteConf.close();
	return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
	ofstream WriteXYZ;

	WriteXYZ.open("frames//config_" + sim_name + "_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "	 " <<  Pbc(y[i]) << "		" << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
		return r - box * rint(r/box);
}

double error(double av, double av2, int n){
	if (n==0) return 0;
	else return sqrt((av2 - av*av)/n);
};

/****************************************************************
*****************************************************************
		_/		_/	_/_/_/	_/			 Numerical Simulation Laboratory
	 _/_/  _/ _/			 _/				Physics Department
	_/	_/_/		_/		_/			 Universita' degli Studi di Milano
 _/		 _/				_/ _/				Prof. D.E. Galli
_/		_/	_/_/_/	_/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/