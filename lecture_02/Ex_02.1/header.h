#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include<stdlib.h>
#include<stdio.h>
#include "random.h"

Random rnd;
int M;                                       // #throws
int N;                                         // #blocks
int block_size;

void Input(void);
double error (double*, double*, int);
