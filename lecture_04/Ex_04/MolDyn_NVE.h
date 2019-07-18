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
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pr;

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
double ave_epot, ave_ekin, ave_pr;
double cum_ave_epot, cum_ave_ekin, cum_ave_pr;
double cum_av2_epot, cum_av2_ekin, cum_av2_pr;
bool restart;
Random rnd;


//functions
void rand_inizialization(void);
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void CorrectOldConf(void);
double Force(int, int);
double Pbc(double);
double error(double,double,int);
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
