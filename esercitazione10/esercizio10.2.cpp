#include <iostream>			// cin, cout: Standard Input/Output Streams Library
#include <fstream>			// Stream class to both read and write from/to files.
#include <cmath>            // Math stuff

#include "mpi.h"

#include "random.h"
#include "parallel.h"

using namespace std;

int main(int argc, char* argv[]){

    // RNG SETUP ##########################################################

    Random *rnd;
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

	// ###################################################################
    // MPI STUFF #########################################################

    MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();

    cout << size <<  " " << rank << endl;

    // ###################################################################
    // CIRCLE WORLD ######################################################

    unsigned const int NCITY = 30;
    unsigned const int ITS  = 60000;
    unsigned const int EVERY = 100;

    const double TMAX = 90;
    const double TMIN = 0;



    double cost, temp;
    vector<unsigned int> best_path;

    vector<vector<double>> radial_cities(NCITY);
    double radius = 30;

    double theta, cityx, cityy;
    for(unsigned int i = 0; i < NCITY; i++){

        theta = rnd->Angle();

        cityx = radius*cos(theta);
        cityy = radius*sin(theta);

        vector<double> city(2);
        
        city[0] = cityx;
        city[1] = cityy;

        radial_cities[i] = city;

    }


    Annealing CircularMap(radial_cities, TMAX, TMIN, ITS);

    cout << "lmao" << endl;


    ofstream fout;
    fout.open("output.circle.dat");
    fout << setw(20) << "Iteration" << setw(20) << "Path Length" << setw(20) << "Temperature" << endl;

    cout << endl << "###################### CIRCULAR MAP ########################" << endl;
    cout << endl << "############################################################" << endl;
    cout << setw(20) << "Iteration" << setw(20) << "Path Length" << setw(20) << "Temperature" << endl;
    for(unsigned int iter = 0; iter < ITS; iter++){

        if (iter % EVERY == 0){

            cost = CircularMap.GetCost();
            temp = CircularMap.LinearScheduler(iter);

            fout << setw(20) << iter << setw(20) << cost << setw(20) << temp << endl;
            //cout << setw(20) << iter << setw(20) << cost << setw(20) << temp << endl;
        }

        CircularMap.Move();

    }
    cout << "############################################################" << endl;
    fout.close();

    fout.open("output.circle.cities.dat");
    fout << setw(20) << "Best Path" << setw(20) << "City X" 
         << setw(20) << "City Y"  << setw(20) << "ID" << endl;

    best_path = CircularMap.GetPath();
    for(unsigned int i = 0; i < radial_cities.size(); i++){

        fout << setw(20) << best_path[i] << setw(20) << radial_cities[i][0] 
             << setw(20) << radial_cities[i][1] << setw(20) << i << endl;

    }
    fout.close();

    cout << endl;

    // ###################################################################
    // SQUARE WORLD ######################################################

    vector<vector<double>> square_cities(NCITY);
    double side = 30;

    for(unsigned int i = 0; i < NCITY; i++){

        cityx = side*rnd->Rannyu();
        cityy = side*rnd->Rannyu();

        vector<double> city(2);
        
        city[0] = cityx;
        city[1] = cityy;

        square_cities[i] = city;

    }

    Annealing SquareMap(square_cities, TMAX, TMIN, ITS);

    fout.open("output.square.dat");
    fout << setw(20) << "Iteration" << setw(20) << "Path Length" << setw(20) << "Temperature" << endl;

    cout << endl << "###################### SQUARE MAP ########################" << endl;
    cout << endl << "############################################################" << endl;
    cout << setw(20) << "Iteration" << setw(20) << "Path Length" << setw(20) << "Temperature" << endl;
    for(unsigned int iter = 0; iter < ITS; iter++){

        if (iter % EVERY == 0){

            cost = SquareMap.GetCost();
            temp = SquareMap.LinearScheduler(iter);

            fout << setw(20) << iter << setw(20) << cost << setw(20) << temp << endl;
            //cout << setw(20) << iter << setw(20) << cost << setw(20) << temp << endl;
        }

        SquareMap.Move();

    }
    cout << "############################################################" << endl;
    fout.close();

    fout.open("output.square.cities.dat");
    fout << setw(20) << "Best Path" << setw(20) << "City X" 
         << setw(20) << "City Y"  << setw(20) << "ID" << endl;

    best_path = SquareMap.GetPath();
    for(unsigned int i = 0; i < radial_cities.size(); i++){

        fout << setw(20) << best_path[i] << setw(20) << square_cities[i][0] 
             << setw(20) << square_cities[i][1] << setw(20) << i << endl;

    }
    fout.close();

    cout << endl;

}