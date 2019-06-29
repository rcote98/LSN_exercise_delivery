
#include "parallel.h"


Annealing :: Annealing(vector< vector<double> >c, double temp_max, double temp_min, unsigned int iterations){

    // RNG SETUP -------------------------------------------------------

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
    
    cities = c;
    ncity = cities.size();
    for(unsigned int i = 0; i < ncity; i++){
        if(!CheckCity(i)){
            cout << "Invalid city vector, exiting..." << endl;
            exit(0);
        }
    }

    iteration = 0;

    path = GeneratePath();

    if(!CheckPath()){
        cout << "Invalid path:" << endl;
        PrintPath();
        cout << "Exiting..." << endl;
        exit(0);
    }


    curr_cost = CostCalc(path);

    tmax = temp_max;
    tmin = temp_min;

    iters = iterations;

}


//###################################################################
// MOVE STUFF #######################################################

void Annealing :: Move(){

    vector<unsigned int> proposed = SwapMutation(path);
    double cost_prop = CostCalc(proposed);

    if(cost_prop < curr_cost){
        path = proposed;
        curr_cost = cost_prop;
    } else {
        
        beta = 1/LinearScheduler(iteration);
        double r = rnd->Rannyu();
        
        if(r < exp(-beta*(cost_prop - curr_cost))){
            path = proposed;
            curr_cost = cost_prop;
        }

    }
    
    iteration++;

}




//###################################################################
// TEMPERATURE SCHEDULERS ###########################################

double Annealing :: LinearScheduler(unsigned int iteration){

    double rate = ((double)iteration/(double)iters);

    return tmax - rate*(tmax-tmin);

}


//###################################################################
// PATH GENERATION ##################################################

vector<unsigned int> Annealing :: GeneratePath() {

    vector<unsigned int> path(ncity);

    for(unsigned int i = 0; i < ncity; i++) path[i] = i;

    for(unsigned int i = 0; i < INIT_MUT; i++){

        path = SwapMutation(path);

    }

    return path;

}

//###################################################################
// SAFETY CHECKS ####################################################

bool Annealing::CheckCity(unsigned int i){ 

    if (cities[i].size() != 2) return false;

    return true;
}

bool Annealing::CheckPath(){
    
    if(path.size() != ncity) return false;

    unsigned int counter = 0;

    for(unsigned int k = 0; k < ncity; k++){
        if (find(path.begin(), path.end(), k) != path.end()) counter++;
    }

    if(counter == ncity) return true;

    return false;
}

//###################################################################
// FITNESS STUFF ####################################################

double Annealing::CostCalc(vector<unsigned int> p){

    double cost = 0;

    for(unsigned int j = 0; j < ncity-1; j++){
        // cost += Norm(p[j], p[j+1]);
        cost += sqrt(Norm(p[j], p[j+1])); 
    }

    // cost += Norm(p[ncity-1], p[0]);
    cost += sqrt(Norm(p[ncity-1], p[0])); 

    return cost;
}


//###################################################################
// MUTATION STUFF ###################################################

 
vector<unsigned int> Annealing::SwapMutation(vector<unsigned int> path){

    int i = (int) (rnd->Rannyu()*ncity);
    int j;

    do{

        j= (int) (rnd->Rannyu()*ncity);

    }while(i == j);

    swap(path[i], path[j]);

    return path;
}



//###################################################################
// AUX FUNCTIONS ####################################################

double Annealing::Norm(unsigned int i, unsigned int j){

    return pow(cities[i][0] - cities[j][0], 2) + pow(cities[i][1] - cities[j][1],2);

}

//###################################################################
// PRINTERS #########################################################


void Annealing::PrintPath(){

    cout << endl << "[";
    for(unsigned int j = 0; j < ncity; j++){
        cout << setw(3) << path[j];
    }
    cout << "]" << endl << endl;


}
