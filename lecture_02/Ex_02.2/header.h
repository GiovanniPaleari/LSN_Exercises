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
int M;                                      // #simulations
int N;                                        // #steps
double step_length=1.;

void Input();
double sq_distance(double, double, double);
double mean(double*, int);
double stddev(double*, int);
int lattice_step(int);
void null(double*, int);
