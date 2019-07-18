#ifndef __GA__h_
#define __GA__h_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "Chromosome.h"
#include "Map.h"

using namespace std;


//Random numbers
#include "random.h"
int seed[4];     //Perchè extern? https://samwho.dev/blog/2013/12/08/duplicate-symbol-what/
Random* rnd;

//Simulation
int npop, ncities, generations;
int SorC;
double l;
double* p_mutation;
double p_crossover;
int mut=5;

Map* cities;
vector<Chromosome> population;
vector<Chromosome> newgen;

ofstream best;

//functions
void Input(void);
void Initialize();
int Selection();
void Crossover(Chromosome, Chromosome,Chromosome*, Chromosome*);
bool order(Chromosome a, Chromosome b);
void Deallocate();


#endif
