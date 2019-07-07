#include <iostream>			// cin, cout: Standard Input/Output Streams Library
#include <fstream>			// Stream class to both read and write from/to files.
#include <cmath>            // Math stuff
#include <vector>
#include <iomanip>
#include <algorithm>

#include "mpi.h"
#include "random.h"

using namespace std;


unsigned const int NCITY = 30;

unsigned const int INIT_MUT = 50;

unsigned const int ITS  = 60000;
unsigned const int EVERY = 100;

const double TMAX = 90;
const double TMIN = 0;

// functions

bool CheckCities(vector<vector<double>> cities);
bool CheckPath(vector<unsigned int> path);
double Cost(vector<unsigned int> path, vector<vector<double>> cities);
vector<unsigned int> SwapMutation(vector<unsigned int> path, Random *rnd);
vector<unsigned int> GeneratePath(Random * rnd);
double LinearScheduler(unsigned int iteration);