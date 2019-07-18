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

bool T;
double step;
int n;

double xtry, ytry, ztry;
double xnow, ynow, znow;
double ave, cum_ave, cum_av2;
double rnow;   //Actual distance
double acc;   //Acceptance probability

void Input(void);
void Metropolis(void);
void Reset(void);
void Measure(void);
void Block(int);
void Equilibration(int);
double psi_1s(double, double, double);
double psi_2p(double, double, double);
void fiftypercent(void);
double error (double, double, int);

#endif
