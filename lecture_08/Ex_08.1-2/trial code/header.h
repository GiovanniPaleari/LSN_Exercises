#ifndef header
#define header

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include "random.h"

using namespace std;

Random rnd;


ifstream input;
ofstream output;

int nstep, nblocks, block_size;
int istep;
int iprint=100;

double step;

double sigma, mu;
double Hpsi;
double xtry;
double xnow;
double ave, cum_ave, cum_av2;
double acc;   //Acceptance probability

void Input(void);
void Metropolis(void);
void Reset(void);
void Measure(void);
void Block(int);
double FTrial(double );
void fiftypercent(void);
/*double sq_distance(double, double, double);
double mean(double*, int);
double stddev(double*, int);
void null(double*, int);*/
double error (double, double, int);

#endif
