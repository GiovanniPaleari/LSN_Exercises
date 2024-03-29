/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "random.h"

//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie, ip, igofr;
double walker[m_props];
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pr;
double err_pot, err_kin, err_etot, err_temp, err_pr, err_gdir;
int nbins;
double bin_size;

// averages
double acc,att;

//configuration
//const int m_part=108;   //cambiare sta roba che m_part != npart
//double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
//double vx[m_part],vy[m_part],vz[m_part];

double *x,*y,*z,*xold,*yold,*zold;
double *vx,*vy,*vz;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblocks, block_size;
int imeasure;
double delta;
bool restart;
Random rnd;
double pi=3.14159265;


//functions
void rand_inizialization(void);
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void CorrectOldConf(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);
void Deallocate(void);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
