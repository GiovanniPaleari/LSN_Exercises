#include "random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <string>

Random rnd;
int M;                                       // #throws
int N;                                         // #blocks
int block_size;
double size;        //size of the horizontal plane
double dist;    //distance between lines
double length;      //lenght of the needle

void Input(void);
double theta ();
double error (double*, double*, int);
