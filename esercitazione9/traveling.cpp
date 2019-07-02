
#include "traveling.h"


TSProblem :: TSProblem(vector< vector<double> > c){

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

    pop_size = ncity*ncity;

    paths.reserve(pop_size);
    new_paths.reserve(pop_size);

    costs.reserve(pop_size);
    fitness.reserve(pop_size);

    generation = 1;

    for(unsigned int i = 0; i < pop_size; i++){
        paths.push_back(GeneratePath());
    }

    for(unsigned int i = 0; i < pop_size; i++){
        if(!CheckPath(i)){
            cout << "Invalid path:" << endl;
            PrintPath(i);
            cout << "Exiting..." << endl;
            exit(0);
        }
    }

    CostCalc();
    SortPopulation();


}


//###################################################################
// GENERATION STUFF #################################################

void TSProblem::AdvanceGeneration(){
  
    unsigned int p1, p2; // parents
    vector<unsigned int> child1, child2; // children
    vector< vector< unsigned int > > children;

    unsigned int moved, shift;
    double rswap, rmswap, rshift, rcross;

    unsigned int npops = 0;

    while(npops < pop_size){

        rcross = rnd->Rannyu();

        if(rcross < PCROSS){
            
            p1 = Selection(SEL_EXP);
            p2 = Selection(SEL_EXP);

            while(p1 == p2){
                p2 = Selection(SEL_EXP);
            }

            children = Crossover(p1,p2);

            child1 = children[0];
            child2 = children[1];


        } else {

            p1 = Selection(SEL_EXP);
            p2 = Selection(SEL_EXP);
            
            while(p1 == p2){
                p2 = Selection(SEL_EXP);
            }

            child1 = paths[p1];
            child2 = paths[p2];

        }

        // swap two cities
        rswap = rnd->Rannyu();
        if (rswap < PSWAP){ 
            child1 = SwapMutation(child1);
            child2 = SwapMutation(child2);
        }

        // swap multiple cities
        rmswap = rnd->Rannyu();
        if (rmswap < PMSWAP){
            moved = (int) (rnd->Rannyu()*ncity/3);
            child1 = SwapMutation(child1, moved);

            moved = (int) (rnd->Rannyu()*ncity/3);
            child2 = SwapMutation(child2, moved);
        }

        // shift the whole path by "shift"
        rshift = rnd->Rannyu();
        if (rshift < PSHIFT){
            shift = (int) (rnd->Rannyu()*ncity/3);
            child1 = ShiftMutation(child1, shift);

            shift = (int) (rnd->Rannyu()*ncity/3);
            child2 = ShiftMutation(child2, shift);
        }

        new_paths[npops] = child1;    
        new_paths[npops+1] = child2;

        //PrintPath(npops, 1);
        //PrintPath(npops+1, 1);

        npops+=2;
    }

    for(unsigned int i = 0; i < pop_size; i++){
        paths[i] = new_paths[i];
    }

    for(unsigned int i = 0; i < pop_size; i++){
        if(!CheckPath(i)){
            cout << "Invalid path: " << i << endl;
            PrintPath(i);
            cout << "Exiting..." << endl;
            exit(0);
        }
    }

    CostCalc();
    SortPopulation();
    generation++;

}


//###################################################################
// PATH GENERATION ##################################################

vector<unsigned int> TSProblem :: GeneratePath() {

    vector<unsigned int> path(ncity);

    for(unsigned int i = 0; i < ncity; i++) path[i] = i;

    for(unsigned int i = 0; i < INIT_MUT; i++){

        path = SwapMutation(path);

    }

    return path;

}

//###################################################################
// SAFETY CHECKS ####################################################

bool TSProblem::CheckCity(unsigned int i){ 

    if (cities[i].size() != 2) return false;

    return true;
}

bool TSProblem::CheckPath(unsigned int i){
    
    if(paths[i].size() != ncity) return false;

    unsigned int counter = 0;

    for(unsigned int k = 0; k < ncity; k++){
        if (find(paths[i].begin(), paths[i].end(), k) != paths[i].end()) counter++;
    }

    if(counter == ncity) return true;

    return false;
}

//###################################################################
// FITNESS STUFF ####################################################

void TSProblem::CostCalc(){

    double cost;
    for(unsigned int i = 0; i < paths.size(); i++){
        cost = 0;

        for(unsigned int j = 0; j < ncity-1; j++){
            //cost += Norm(paths[i][j], paths[i][j+1]);
            cost += sqrt(Norm(paths[i][j], paths[i][j+1]));
        }
        cost += sqrt(Norm(paths[i][ncity-1], paths[i][0]));
        
        costs[i] = cost;
        fitness[i] = 1/cost;
    }
}

void TSProblem::SortPopulation(){
    QuickSort(0, pop_size-1);
}

//###################################################################
// MUTATION STUFF ###################################################

unsigned int TSProblem::Selection(double exponent){

    return (int)(pop_size*pow(rnd->Rannyu(),exponent));

}

vector< vector< unsigned int > > TSProblem::Crossover(unsigned int p1, unsigned int p2){

    unsigned int cut = (int)(ncity*rnd->Rannyu());
    unsigned int tail = ncity - cut;

    vector<unsigned int> child1(ncity), child2(ncity);
    vector<unsigned int> missing1(tail), missing2(tail);

    unsigned int count;

    for(unsigned int i = 0; i<cut; i++){
        child1[i] = paths[p1][i];
        child2[i] = paths[p2][i];
    }

    for(unsigned int i = 0; i < tail; i++){
        missing1[i] = paths[p1][cut + i];
        missing2[i] = paths[p2][cut + i];
    }

    count = 0;
    while(count < tail){
        for(unsigned int i=0; i< ncity; i++){
            for(unsigned int j=0; j<tail; j++){
	            if(paths[p2][i]==missing1[j]){
	                child1[cut+count]=missing1[j];
	                count++;
	            }
            }
        }
    }

    count = 0;
    while(count < tail){
        for(unsigned int i=0; i< ncity; i++){
            for(unsigned int j=0; j<tail; j++){
	            if(paths[p1][i]==missing2[j]){
	                child2[cut+count]=missing2[j];
	                count++;
	            }
            }
        }
    }

    vector< vector< unsigned int > > children;

    children.push_back(child1);
    children.push_back(child2);

    return children;
}

vector<unsigned int> TSProblem::SwapMutation(vector<unsigned int> path){

    int i = (int) (rnd->Rannyu()*ncity);
    int j;

    do{

        j= (int) (rnd->Rannyu()*ncity);

    }while(i == j);

    swap(path[i], path[j]);

    return path;
}

vector<unsigned int> TSProblem::SwapMutation(vector<unsigned int> path, unsigned int m){

    int i = (int) (rnd->Rannyu()*(ncity - m));
    int j;

    do{

        j=(int)(rnd->Rannyu()*(ncity - m));

    }while(i == j && abs(j-i) <= m);

    for(unsigned int k = 0; k < m; k++){
        swap(path[i+k], path[j+k]);
    }
    

    return path;
}

vector<unsigned int> TSProblem::ShiftMutation(vector<unsigned int> path, unsigned int m){

    double cache[m];

    for(unsigned int i = 0; i < m; i++){
        cache[i] = path[i];
    }
    for(unsigned int i = 0; i < ncity-m; i++){
        path[i] = path[i+m];
    }
    for(unsigned int i = 0; i < m; i++){
        path[ncity-m+i] = cache[i];
    }

    return path;
}

//vector<unsigned int> TSProblem::MShiftMutation(vector<unsigned int> path, unsigned int m, unsigned int n){


  //  return void;
//}


//###################################################################
// STATISTICS STUFF #################################################

double TSProblem::AverageFitness(){

    double sum = 0;
    for(unsigned int i = 0; i < pop_size; i++){
        sum += fitness[i];
    }
    return sum/pop_size;

}


double TSProblem::HalfAverageFitness(){

    double sum = 0;
    for(unsigned int i = 0; i < pop_size/2; i++){
        sum += fitness[i];
    }
    return sum/(pop_size/2);

}

double TSProblem::AveragePathLength(){

    double sum = 0;
    for(unsigned int i = 0; i < pop_size; i++){
        sum += costs[i];
    }
    return sum/pop_size;
}

double TSProblem::HalfAveragePathLength(){

    double sum = 0;
    for(unsigned int i = 0; i < pop_size/2; i++){
        sum += costs[i];
    }
    return sum/(pop_size/2);

}



//###################################################################
// AUX FUNCTIONS ####################################################

double TSProblem::Norm(unsigned int i, unsigned int j){

    return pow(cities[i][0] - cities[j][0], 2) + pow(cities[i][1] - cities[j][1],2);

}

int TSProblem::Partition(int low, int high) 
{ 
    int pivot = costs[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than or 
        // equal to pivot 
        if (costs[j] <= pivot) 
        { 
            i++;    // increment index of smaller element 

            swap(costs[i], costs[j]); 
            swap(paths[i], paths[j]);
            swap(fitness[i], fitness[j]); 
        } 
    } 

    swap(costs[i + 1], costs[high]); 
    swap(paths[i + 1], paths[high]);
    swap(fitness[i+1], fitness[high]); 
    return (i + 1); 
} 
  
void TSProblem::QuickSort(int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = Partition(low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        QuickSort(low, pi - 1); 
        QuickSort(pi + 1, high); 
    } 
} 
//###################################################################
// PRINTERS #########################################################

void TSProblem::ShowPops(unsigned int n){

    cout << endl;
    cout << setw(15) << "GENERATION: " << setw(15) << generation << endl;
    cout << setw(15) << "SAMPLE SIZE: " << setw(15) << n << endl;
    cout << setw(15) << "POP SIZE: " << setw(15) << pop_size << endl;
    cout << setw(15) << "AVG FIT: " << setw(15) << AverageFitness() << endl; 
    cout << setw(15) << "AVG DIST: " << setw(15) << AveragePathLength() << endl;
    cout << endl;

    cout << setw(92) << "PATH" << setw(15) << "DISTANCE" << setw(15) << "FITNESS" << endl;
    for(unsigned int i = 0; i < n; i++){
        cout << "[";
        for(unsigned int j = 0; j < ncity; j++){
            cout << setw(3) << paths[i][j];
        }
        cout << "]" << setw(15) << costs[i] << setw(15) << fitness[i] << endl;
    }

}

void TSProblem::PrintPath(unsigned int i, unsigned int newpaths){

    if(newpaths == 0){
        cout << endl << "[";
        for(unsigned int j = 0; j < ncity; j++){
            cout << setw(3) << paths[i][j];
        }
        cout << "]" << endl << endl;
    } else {
        cout << endl << "[";
        for(unsigned int j = 0; j < ncity; j++){
            cout << setw(3) << new_paths[i][j];
        }
        cout << "]" << endl << endl;
    }

}
