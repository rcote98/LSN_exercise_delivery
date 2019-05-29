#include <iostream>			// cin, cout: Standard Input/Output Streams Library
#include <fstream>			// Stream class to both read and write from/to files.
#include <cmath>            // Math stuff

#include "random.h"
#include "traveling.h"

using namespace std;

int main(){

    // RNG SETUP -------------------------------------------------------

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

	// ------------------------------------------------------------------

    unsigned const int NCITY = 30;
    unsigned const int GENS = 100;

    int gen;
    double avg_path;
    double ave_fit;
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


    TSProblem CircularMap(radial_cities);

    ofstream fout;
    fout.open("output.circle.gen.dat");
    fout << setw(20) << "Generation" << setw(20) << "Avg. Path" << setw(20) << "Avg. Fitness" << endl;

    cout << endl << "###################### CIRCULAR MAP ########################" << endl;
    cout << endl << "############################################################" << endl;
    cout << setw(20) << "Generation" << setw(20) << "Avg. Path" << setw(20) << "Avg. Fitness" << endl;
    for(unsigned int i = 0; i < GENS; i++){

        gen = CircularMap.GetGeneration();
        avg_path = CircularMap.HalfAveragePathLength();
        ave_fit = CircularMap.HalfAverageFitness();

        fout << setw(20) << gen << setw(20) << avg_path << setw(20) << ave_fit << endl;

        cout << setw(20) << gen << setw(20) << avg_path << setw(20) << ave_fit << endl;

        CircularMap.AdvanceGeneration();

    }
    cout << "############################################################" << endl;
    fout.close();

    fout.open("output.circle.cities.dat");
    fout << setw(20) << "City ID" << setw(20) << "City X" 
         << setw(20) << "City Y"  << setw(20) << "Best Path" << endl;
    best_path = CircularMap.GetBestPath();
    for(unsigned int i = 0; i < radial_cities.size(); i++){

        fout << setw(20) << i << setw(20) << radial_cities[i][0] 
             << setw(20) << radial_cities[i][1] << setw(20) << best_path[i] << endl;

    }
    fout.close();

    CircularMap.ShowPops(10);

    cout << endl;

    //#########################################################################

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


    TSProblem SquareMap(square_cities);

    fout.open("output.square.gen.dat");
    fout << setw(20) << "Generation" << setw(20) << "Avg. Path" << setw(20) << "Avg. Fitness" << endl;

    cout << endl << "###################### SQUARE MAP ########################" << endl;
    cout << endl << "############################################################" << endl;
    cout << setw(20) << "Generation" << setw(20) << "Avg. Path" << setw(20) << "Avg. Fitness" << endl;
    for(unsigned int i = 0; i < GENS; i++){

        gen = SquareMap.GetGeneration();
        avg_path = SquareMap.HalfAveragePathLength();
        ave_fit = SquareMap.HalfAverageFitness();

        fout << setw(20) << gen << setw(20) << avg_path << setw(20) << ave_fit << endl;

        cout << setw(20) << gen << setw(20) << avg_path << setw(20) << ave_fit << endl;

        SquareMap.AdvanceGeneration();

    }
    cout << "############################################################" << endl;
    fout.close();

    fout.open("output.square.cities.dat");
    fout << setw(20) << "City ID" << setw(20) << "City X" 
         << setw(20) << "City Y"  << setw(20) << "Best Path" << endl;
    best_path = SquareMap.GetBestPath();
    for(unsigned int i = 0; i < square_cities.size(); i++){

        fout << setw(20) << i << setw(20) << square_cities[i][0] 
             << setw(20) << square_cities[i][1] << setw(20) << best_path[i] << endl;

    }
    fout.close();

    SquareMap.ShowPops(10);

    cout << endl;
}