#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include "random.h"

Random rnd;
int N;
double decay_rate=1.;
double width=1., median=0;
int nbin=100;
int h_std[100], h_exp[100], h_lrz[100];		//vectors representing histograms

void Input(void);
void ResetHisto(void);
void fill_histo (int nbin, double min, double max, int* h, double x);
