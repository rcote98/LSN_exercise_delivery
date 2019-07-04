/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "esercizio4.h"

using namespace std;


int main(int argc, char * argv[]){ 

  verbose = false;

  arg_count = argc;
	args = argv;

  Input();             //Inizialization
  int nconf = 1;
  for(unsigned int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%measure_step == 0){
        Measure();     //Properties measurement
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  Averages();          //Data blocking implementation
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation

	// RNG SETUP ----------------------------------------------

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

	// ----------------------------------------------------------
  // Simulation Setup -----------------------------------------

  ifstream ReadInput,ReadConf;

  if(verbose){ // you know the drill...
    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "The program uses Lennard-Jones units " << endl;
  }

  //Read Input File
  ReadInput.open("input.dat"); 
  ReadInput >> temp;
  ReadInput >> npart;
  ReadInput >> rho;
  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblocks;
  ReadInput >> iprint;
  ReadInput >> RESTART;
  ReadInput.close();

  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  max_radius = box/2;


  if (arg_count == 2){

		temp = atof(args[1]);

	}

  // Simulation Config
  cout << "Temperature = " << temp << endl;
  cout << "Number of particles = " << npart << endl;
  cout << "Density of particles = " << rho << endl;
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  if(RESTART) cout << "Starting simulation from zero..." << endl << endl;
  else cout << "Resuming simulation from previous state..." << endl << endl; 
  

  if (RESTART) { // Simulation beginning from zero

    // Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (unsigned int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

    //Prepare initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};

    for (unsigned int i=0; i<npart; ++i){
      vx[i] = rnd->Rannyu() - 0.5;
      vy[i] = rnd->Rannyu() - 0.5;
      vz[i] = rnd->Rannyu() - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }

    for (unsigned int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

    double sumv2 = 0.0, fs;
    for (unsigned int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (unsigned int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = x[i] - vx[i] * delta;
      yold[i] = y[i] - vy[i] * delta;
      zold[i] = z[i] - vz[i] * delta;
    }


  } else { // Simulation continuing from restored state

    // Read x(t) configuration
    cout << "Reading x(t) configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (unsigned int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

    // Read x(t-dt) configuration
    cout << "Reading x(t-dt) configuration from file config.0.prev " << endl << endl;
    ReadConf.open("config.0.prev");
    for (unsigned int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf.close();

    // Calculate x(t+dt) and v(t)


    double xnew[npart], ynew[npart], znew[npart]; 
    double vxnew[npart], vynew[npart], vznew[npart]; 
    double fx[m_part], fy[m_part], fz[m_part];

    unsigned int recalcs = 500; // recalcs of the previous position

    for(unsigned int k = 0; k < recalcs; k++){

      for(unsigned int i=0; i<npart; ++i){ //Force acting on particle i
        fx[i] = Force(i,0);
        fy[i] = Force(i,1);
        fz[i] = Force(i,2);
      }

      for(unsigned int i=0; i<npart; ++i){ //Verlet integration scheme

        xnew[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

        vxnew[i] = Pbc(xnew[i] - xold[i])/(2.0 * delta);
        vynew[i] = Pbc(ynew[i] - yold[i])/(2.0 * delta);
        vznew[i] = Pbc(znew[i] - zold[i])/(2.0 * delta);
      }

      double sumv[3] = {0.0, 0.0, 0.0};
      for (unsigned int i=0; i<npart; ++i){
        sumv[0] += vxnew[i];
        sumv[1] += vynew[i];
        sumv[2] += vznew[i];
      }
      
      for (unsigned int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

      double sumv2 = 0.0, fs;
      for (unsigned int i=0; i<npart; ++i){
        vxnew[i] = vxnew[i] - sumv[0];
        vynew[i] = vynew[i] - sumv[1];
        vznew[i] = vznew[i] - sumv[2];

        sumv2 += vxnew[i]*vxnew[i] + vynew[i]*vynew[i] + vznew[i]*vznew[i];
      }

      sumv2 /= (double)npart;

      fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
      for (unsigned int i=0; i<npart; ++i){
        vxnew[i] *= fs;
        vxnew[i] *= fs;
        vznew[i] *= fs;

        xold[i] = xnew[i] - vxnew[i]*2*delta;
        yold[i] = ynew[i] - vynew[i]*2*delta;
        zold[i] = znew[i] - vznew[i]*2*delta;
      }

    }

    Move();

    double k=0.0;
    for (unsigned int i=0; i<npart; ++i) k += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    stima_temp = (2.0 / 3.0) * k/(double)npart;

    cout << "Temperature after adjustment: " << stima_temp << endl << endl;

  }
   return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew; 
  double fx[m_part], fy[m_part], fz[m_part];

  for(unsigned int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(unsigned int i=0; i<npart; ++i){ //Verlet integration scheme

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

double Force(unsigned int ip, unsigned int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (unsigned int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
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
  double v, p, t, r, vij, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres, Gofr;

  Epot.open("output.epot.dat",ios::app);
  Ekin.open("output.ekin.dat",ios::app);
  Temp.open("output.temp.dat",ios::app);
  Etot.open("output.etot.dat",ios::app);
  Pres.open("output.pres.dat",ios::app);
  Gofr.open("output.gofr.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0; 
  p = 0.0;
  r = 0.0;

  //cycle over pairs of particles
  for (unsigned int i=0; i<npart-1; ++i){
    for (unsigned int j=i+1; j<npart; ++j){

      dx = Pbc( x[i] - x[j] );
      dy = Pbc( y[i] - y[j] );
      dz = Pbc( z[i] - z[j] );

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      

      // g(r)
      for(unsigned int k = 0; k < nbins; k++){
          interval_min = k*max_radius/nbins;
          interval_max = (k+1)*max_radius/nbins;

          norm = rho*npart*4*M_PI/3*(pow(interval_max, 3) - pow(interval_min, 3));

          if(dr < interval_max && dr > interval_min){
            gofr_hist[k] += 2./norm;
          }
      }

      if(dr < rcut){
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        pij = 48.0*(pow(dr,-12) - 0.5*pow(dr,-6));

        // Potential energy
        v += vij;

        // Pressure
        p += pij;

      }
    }
  }  


  //Kinetic energy
  for (unsigned int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    
  stima_pot = v/(double)npart; //Potential energy
  stima_kin = t/(double)npart; //Kinetic energy
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy
  stima_pres = (rho*(2.0/3.0)*t + p/(3*v))/(double)npart; // Pressure

	stima_gofr = 0; // g(r) average
	for(unsigned int k = 0; k < nbins; k++){

		interval_min = k*max_radius/nbins;
		interval_max = (k+1)*max_radius/nbins;
		
		r = (interval_max - interval_min)/2;

		stima_gofr += r*gofr_hist[k];

    // reset histogram
    gofr_hist[k] = 0;

	}

  

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  Pres << stima_pres << endl;
  Gofr << stima_gofr << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Pres.close();
  Gofr.close();

  return;

}



void Averages(void){


	ifstream Epot, Ekin, Etot, Temp, Pres, Gofr;
  double *etot_ave, *epot_ave, *ekin_ave;
  double *temp_ave, *pres_ave, *gofr_ave; 

	Epot.open("output.epot.dat");
	Ekin.open("output.ekin.dat");
	Temp.open("output.temp.dat");
	Etot.open("output.etot.dat");
	Pres.open("output.pres.dat");
  Gofr.open("output.gofr.dat");

	unsigned int nmeasures = nstep/measure_step;

	unsigned int k;
	double L = nmeasures/nblocks;

  unsigned int nprops = 6;
	double* sum = new double[nprops];


	double etot_meas[nmeasures];
	double ekin_meas[nmeasures];
	double epot_meas[nmeasures];
	double temp_meas[nmeasures];
	double pres_meas[nmeasures];
  double gofr_meas[nmeasures];

	for(unsigned int i = 0; i < nmeasures; i++)
	{
			Epot >> epot_meas[i];
			Ekin >> ekin_meas[i];
			Temp >> temp_meas[i];
			Etot >> etot_meas[i];
			Pres >> pres_meas[i];
      Gofr >> gofr_meas[i];
	}
	
	etot_ave = new double[nblocks];
	epot_ave = new double[nblocks];
	ekin_ave = new double[nblocks];
	temp_ave = new double[nblocks];
	pres_ave = new double[nblocks];
  gofr_ave = new double[nblocks];

	for(unsigned int i = 0; i<nblocks; i++){
		
		for(unsigned int j = 0; j < nprops; j++) sum[j] = 0;


		for(unsigned int j = 0; j < L; j++){

			k = j+i*L;

			sum[0] += epot_meas[k];
			sum[1] += ekin_meas[k];
			sum[2] += temp_meas[k];
			sum[3] += etot_meas[k];
			sum[4] += pres_meas[k];
      sum[5] += gofr_meas[k];
		}

		epot_ave[i] = sum[0]/L;
		ekin_ave[i] = sum[1]/L;
		temp_ave[i] = sum[2]/L;
		etot_ave[i] = sum[3]/L;
		pres_ave[i] = sum[4]/L;
    gofr_ave[i] = sum[5]/L;


	}

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();
	Pres.close();
  Gofr.close();

	data_blocking(nblocks, etot_ave, "output.etot_ave.dat");
	data_blocking(nblocks, epot_ave, "output.epot_ave.dat");
	data_blocking(nblocks, ekin_ave, "output.ekin_ave.dat");
	data_blocking(nblocks, temp_ave, "output.temp_ave.dat");
	data_blocking(nblocks, pres_ave, "output.pres_ave.dat");
  data_blocking(nblocks, gofr_ave, "output.gofr_ave.dat");

  cout << endl << "Final temperature " << temp_ave[nblocks-1] << endl << endl;

  stima_temp = temp_ave[nblocks-1];

}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Printing x(t-dt) configuration to file config.final.prev " << endl << endl;
  WriteConf.open("config.final.prev");

  for (unsigned int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();

  cout << "Printing x(t) configuration to file config.final " << endl;
  WriteConf.open("config.final");

  for (unsigned int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  cout << "Printing final temperature to temp.final" << endl;
	WriteConf.open("temp.final");

	WriteConf << stima_temp << endl;
	
	WriteConf.close();


  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (unsigned int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void data_blocking(unsigned int N, double* ave, string fname){

	double * sum_prog = new double[N];
	double * su2_prog = new double[N];
	double * err_prog = new double[N];

	double av;
	for(unsigned int i = 0; i<N; i++){
		for(unsigned int j=0; j<i+1; j++){
			av = ave[j];
			sum_prog[i] += av;
			su2_prog[i] += pow(av,2);
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		su2_prog[i] = su2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog[i], su2_prog[i] ,i);
	}

	fstream fout;
	fout.open(fname, fstream::out);
	fout << "block, sum_prog, err_prog" << endl;
	for(unsigned int i = 0; i < N; i++){
		fout << i+1 << ", " << sum_prog[i] << ", " << err_prog[i] << endl;
	}
	fout.close();

	delete sum_prog;
	delete su2_prog;
	delete err_prog;

}

double error(double av, double av2, int n){
	if (n==0) return 0;
	else return sqrt((av2 - av*av)/n);
};

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/