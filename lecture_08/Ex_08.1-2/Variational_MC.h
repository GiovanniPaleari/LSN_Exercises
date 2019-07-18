/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __var_
#define __var_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;


// averages
double Hpsi;
double blk_av, blk_norm,accepted,attempted;
double glob_av, glob_av2;
double stima, err;

// simulation
int nstep, nblk, mcstep;
double jump;
double sigma, mu;
double xin, xnow, xtry;

int nbin = 100;
double bin_size;
double histo[100];

//pigreco
const double pi=3.141592654;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
double Error(double,double,int);
double FTrial(double);
void fiftypercent(void);
void PrintHisto(void);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
