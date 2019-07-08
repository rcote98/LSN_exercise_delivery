#include "esercizio10.2.h"

using namespace std;

int main(int argc, char* argv[]){

    // ###################################################################
    // MPI STUFF #########################################################

    MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();

    // ####################################################################
    // RNG SETUP ##########################################################

    Random *rnd;
	rnd = new Random();

    int seed[4];
    int * ps = new int[size * 2];
	int * sub_ps = new int[2];

	if(rank == 0){

            // Manager reads the primes and assigns them to slaves

		    ifstream Primes("Primes");

		    string primes_input;
		
            if(Primes.is_open())
            {
                // Assign to each subprocess a different line of the file,
                // so that the generated numbers are different.
                for(int i = 0; i < size; ++i)
                {
                    Primes >> ps[2*i] >> ps[2*i + 1];
                }
            }
            else{
                cerr << "Unable to open Primes." << std::endl;
                Primes.close();
            }

	    }

        MPI_Scatter(
                ps,             // Array of data to be sent
                2,              // Length of sent data
                MPI_INT,        // Type of sent data
                sub_ps,         // Address to hold received data
                2,              // Length of received data
                MPI_INT,        // Type of received data
                0,              // Rank of manager process
                MPI::COMM_WORLD // Communicator
        );


    ifstream input("seed.in");
    string property;
    if(input.is_open())
    {
        if(!input.eof())
        {
            input >> property;
            if(property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd->SetRandom(seed, sub_ps[0], sub_ps[1]);
            }
        }
        input.close();
    }
    else
    {
        cerr << "Unable to open seed.in." << std::endl;
    }


    // ###################################################################
    // CITY CREATION #####################################################

    vector<vector<double>> radial_cities(NCITY);
    double radius = 30;
    vector<vector<double>> square_cities(NCITY);
    double side = 30;


    if(rank == 0){

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

        for(unsigned int i = 0; i < NCITY; i++){

            cityx = side*rnd->Rannyu();
            cityy = side*rnd->Rannyu();

            vector<double> city(2);
            
            city[0] = cityx;
            city[1] = cityy;

            square_cities[i] = city;

        }



    }

    // ###################################################################
    // CITY BROADCAST ####################################################


    double circle_cities_x[NCITY];
    double circle_cities_y[NCITY];

    double square_cities_x[NCITY];
    double square_cities_y[NCITY];

    if(rank == 0)
	{
		for(unsigned int i = 0; i < NCITY; ++i){

			circle_cities_x[i] = radial_cities[i][0];
			circle_cities_y[i] = radial_cities[i][1];

            square_cities_x[i] = square_cities[i][0];
			square_cities_y[i] = square_cities[i][1];
		}
	}


	MPI_Bcast(
			circle_cities_x,        // Pointer to the beginning of the array
			NCITY,                  // Length of the array
			MPI_DOUBLE,             // Type of the array values
			0,                      // Manager process rank
			MPI::COMM_WORLD         // Communicator
		);

	MPI_Bcast(
			circle_cities_y,        // Pointer to the beginning of the array
			NCITY,                  // Length of the array
			MPI_DOUBLE,             // Type of the array values
			0,                      // Manager process rank
			MPI::COMM_WORLD         // Communicator
        );

    MPI_Bcast(
			square_cities_x,        // Pointer to the beginning of the array
			NCITY,                  // Length of the array
			MPI_DOUBLE,             // Type of the array values
			0,                      // Manager process rank
			MPI::COMM_WORLD         // Communicator
		);

	MPI_Bcast(
			square_cities_y,        // Pointer to the beginning of the array
			NCITY,                  // Length of the array
			MPI_DOUBLE,             // Type of the array values
			0,                      // Manager process rank
			MPI::COMM_WORLD         // Communicator
        );



    if(rank > 0)
	{
		// Unpack into the nice format again

        for(unsigned int i = 0; i < NCITY; i++){

            vector<double> empty = {0,0};

            radial_cities[i] = empty;
			radial_cities[i][0] = circle_cities_x[i];
            radial_cities[i][1] = circle_cities_y[i];

            square_cities[i] = empty;
            square_cities[i][0] = square_cities_x[i];
            square_cities[i][1] = square_cities_y[i];

		}
    }

    MPI_Barrier(MPI::COMM_WORLD);

    // ###################################################################
    // PATH GENERATION ###################################################

    
    vector<unsigned int> circle_path(NCITY);
    vector<unsigned int> square_path(NCITY);

    circle_path = GeneratePath(rnd);
    square_path = GeneratePath(rnd);

    double circle_cost = Cost(circle_path, radial_cities);
    double square_cost = Cost(square_path, square_cities);

    // ###################################################################
    // SIMULATED ANNEALING ###############################################

    vector<unsigned int> proposed_c;
    vector<unsigned int> proposed_s;
    double cost_prop_c, cost_prop_s;

    double temp, r;

    ofstream foutc, fouts;
    foutc.open("output.circle.p" +  to_string(rank) + ".dat", fstream::out);
    fouts.open("output.square.p" +  to_string(rank) + ".dat", fstream::out);

    foutc << setw(20) << "Iteration" << setw(20) << "Path Length" << setw(20) << "Temperature" << endl;
    fouts << setw(20) << "Iteration" << setw(20) << "Path Length" << setw(20) << "Temperature" << endl;

    for(unsigned int iter = 0; iter < ITS; iter++){

        temp = LinearScheduler(iter);

        proposed_c = SwapMutation(circle_path, rnd);
        cost_prop_c = Cost(proposed_c, radial_cities);

        if(cost_prop_c < circle_cost){

            circle_path = proposed_c;
            circle_cost = cost_prop_c;

        } else {
        
            r = rnd->Rannyu();    
            if(r < exp(-(cost_prop_c - circle_cost)/temp)){

                circle_path = proposed_c;
                circle_cost = cost_prop_c;
            }

        }

        proposed_s = SwapMutation(square_path, rnd);
        cost_prop_s = Cost(proposed_s, square_cities);

        if(cost_prop_s < square_cost){

            square_path = proposed_s;
            square_cost = cost_prop_s;

        } else {
        
            r = rnd->Rannyu();
            if(r < exp(-(cost_prop_s - square_cost)/temp)){

                square_path = proposed_s;
                square_cost = cost_prop_s;
            }
            
        }

        if (iter % EVERY == 0){

            foutc << setw(20) << iter << setw(20) << circle_cost << setw(20) << temp << endl;
            fouts << setw(20) << iter << setw(20) << square_cost << setw(20) << temp << endl;

        }

    }

    foutc.close();
    fouts.close();

    cout << "(" << rank << "): Annealing Finished!"  << endl;

    MPI_Barrier(MPI::COMM_WORLD);

    double best_circle[2], best_square[2]; 
 
    MPI_Allreduce(
				&circle_cost,
				&best_circle,
				1,
				MPI_LONG_DOUBLE_INT,
				MPI_MINLOC,
				MPI::COMM_WORLD
    );

    MPI_Allreduce(
				&square_cost,
				&best_square,
				1,
				MPI_LONG_DOUBLE_INT,
				MPI_MINLOC,
				MPI::COMM_WORLD
    );

	if(rank == 3){

        foutc.open("output.bestcircle.cities.dat");
        foutc << setw(20) << "Best Path" << setw(20) << "City X" 
              << setw(20) << "City Y"  << setw(20) << "ID" << endl;

        for(unsigned int i = 0; i < NCITY; i++){
            foutc << setw(20) << circle_path[i] << setw(20) << radial_cities[i][0] 
             << setw(20) << radial_cities[i][1] << setw(20) << i << endl;
        }

        foutc.close();
    
    }

    if(rank == 1){

        fouts.open("output.bestsquare.cities.dat");
        fouts << setw(20) << "Best Path" << setw(20) << "City X" 
              << setw(20) << "City Y"  << setw(20) << "ID" << endl;

        for(unsigned int i = 0; i < NCITY; i++){
            fouts << setw(20) << square_path[i] << setw(20) << square_cities[i][0] 
             << setw(20) << square_cities[i][1] << setw(20) << i << endl;
        }

        fouts.close();

    }

    MPI::Finalize();

};

// ###################################################################
// FUNCTION DEFINITION ###############################################

bool CheckCities(vector<vector<double>> cities){ 

    for(unsigned int i = 0; i < cities.size(); i++){
        if (cities[i].size() != 2) return false;
        if (cities[i][0] == 0 && cities[i][1] == 0) return false;
    }

    return true;
}

bool CheckPath(vector<unsigned int> path){
    
    if(path.size() != NCITY) return false;

    unsigned int counter = 0;

    for(unsigned int k = 0; k < NCITY; k++){
        if (find(path.begin(), path.end(), k) != path.end()) counter++;
    }

    if(counter == NCITY) return true;

    return false;
}


double Cost(vector<unsigned int> path, vector<vector<double>> cities){

    double cost = 0;

    for(unsigned int j = 0; j < NCITY-1; j++){
        cost += sqrt( pow(cities[path[j]][0] - cities[path[j+1]][0], 2) + pow(cities[path[j]][1] - cities[path[j+1]][1], 2)); 
    }
    cost += sqrt(pow(cities[path[0]][0] - cities[path[NCITY-1]][0], 2) + pow(cities[path[0]][1] - cities[path[NCITY-1]][1], 2)); 

    return cost;
}

vector<unsigned int> SwapMutation(vector<unsigned int> path, Random *rnd){

    int a = (int) (rnd->Rannyu()*NCITY);
    int b;

    do{

        b= (int) (rnd->Rannyu()*NCITY);

    }while(a == b);

    swap(path[a], path[b]);

    return path;
}

vector<unsigned int> GeneratePath(Random * rnd){

    vector<unsigned int> path(NCITY);

    for(unsigned int i = 0; i < NCITY; i++) path[i] = i;

    for(unsigned int i = 0; i < INIT_MUT; i++){
        path = SwapMutation(path, rnd);
    }

    return path;
}

double LinearScheduler(unsigned int iteration){

    double rate = ((double)iteration/(double)ITS);

    return TMAX - rate*(TMAX-TMIN);

}