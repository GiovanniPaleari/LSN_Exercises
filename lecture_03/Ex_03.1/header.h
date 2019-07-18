#ifndef header
#define header

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include "random.h"

using namespace std;

Random rnd;

FILE* out;
ifstream input;
ofstream output;

int M=10000;                                      // #simulations
int N=100;                                        // #blocks
int block_size=M/N;

double ave[2];
double cum_ave[2], cum_av2[2];
int icall=0, iput=1;
int S_0;																			//asset price at t=0
double T;																			//delivery time
int K;																				//strike price
double r;																			//risk free interest rate
double sigma;																  //volatility


void Reset(void);
void Input(void);
double error (double*, double*, int);
void call(int,int,int,string);
void put(int,int,int,string);

#endif
