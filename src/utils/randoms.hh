#ifndef RANDOMS_UTIL_HH
#define RANDOMS_UTIL_HH


#include <random>
#include <ctime>
#include <cstdlib>

int randint (int i) { return std::rand()%i;}
int randint (int i,int j) { return (std::rand()%(j-i))+i;}


#endif