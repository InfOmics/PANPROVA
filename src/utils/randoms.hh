#ifndef RANDOMS_UTIL_HH
#define RANDOMS_UTIL_HH


#include <random>
#include <ctime>
#include <cstdlib>      // std::rand, std::srand
#include <algorithm>


// random number generator
std::mt19937_64 rng;
std::uniform_real_distribution<double> unif(0, 1);

// initialize the random seed for the random number generator
void init_random_seed(int seed){
    std::srand ( unsigned ( seed ) );
    rng.seed(seed);
}

// returns a random integer in the range [0, i)
int randint (int i) { return std::rand()%i;}

// returns a random integer in the range [i, j)
int randint (int i,int j) { return (std::rand()%(j-i))+i;}

double generate_unif(){
    return unif(rng);
}

#endif