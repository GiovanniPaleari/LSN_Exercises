/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  Save();

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  mc_time=0;

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

  if (restart == 0){
//initial configuration: T=infinite
    for (int i=0; i<nspin; ++i){
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }else{
    ReadConf.open("config.final");
    for (int i=0; i<nspin; ++i){
      ReadConf >> s[i];
    }
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro){
  int o;
  double p, energy_old, energy_new;
//  double energy_up, energy_down, sm;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    energy_old = Boltzmann(s[o],o);
    energy_new = Boltzmann(-1*s[o],o);

    if(metro==1) //Metropolis
    {
        p = exp(-beta*(energy_new-energy_old));
        double r = rnd.Rannyu();
        if (r<p){s[o] *= -1; accepted++;}
    }
    else //Gibbs sampling
    {
      p=1./(double)(1+exp(beta*(energy_new-energy_old)));    //compute conditional probability
      double r = rnd.Rannyu();
      if(r<p){s[o] *= -1; accepted++;}
    }
    attempted++;
  }
  mc_time++;

}

double Boltzmann(int sm, int ip)      //sm = current spin
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;      //if we switch spin sm, this is the only part of energy which changes
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }

  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;

}


void Reset(int iblk) //Reset block averages
{

  if(iblk == 1){
    for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
      }
  }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    string energy = "output.ene_";
    string spheat = "output.heat_";
    string magn = "output.mag_";
    string susc = "output.chi_";
    string ext = ".0";
    string t = to_string(temp);
    t.erase (3, std::string::npos);
    //t.erase (t.find_last_not_of('0') + 1, std::string::npos);
    string sampling;

    if (metro==1){ sampling="Metro_";}else{sampling="Gibbs_";}

    if (h==0){

      Ene.open(sampling + energy + t + ext,ios::app);
      stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
      glob_av[iu]  += stima_u;
      glob_av2[iu] += stima_u*stima_u;
      err_u=Error(glob_av[iu],glob_av2[iu],iblk);
      Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.close();

      Heat.open(sampling + spheat + t + ext,ios::app);
      stima_c = pow(beta,2.)*(blk_av[ic]/blk_norm - pow(stima_u*nspin, 2.))/(double)nspin; //Specific Heat
      glob_av[ic]  += stima_c;
      glob_av2[ic] += stima_c*stima_c;
      err_c=Error(glob_av[ic],glob_av2[ic],iblk);
      Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      Heat.close();

      Chi.open(sampling + susc + t + ext,ios::app);
      stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
      glob_av[ix]  += stima_x;
      glob_av2[ix] += stima_x*stima_x;
      err_x=Error(glob_av[ix],glob_av2[ix],iblk);
      Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();

    }else{

      Mag.open(sampling + magn + "h0.02_T" + t + ext,ios::app);
      stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
      glob_av[im]  += stima_m;
      glob_av2[im] += stima_m*stima_m;
      err_m=Error(glob_av[im],glob_av2[im],iblk);
      Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();
    }

    cout << "----------------------------" << endl << endl;
}

void Save(void){

  string name1 = "UCXofT.txt";
  string name2 = "MofT.txt";
  string sampling;

  if (metro==1){ sampling="Metro_";}else{sampling="Gibbs_";}

  if (restart==1){
    if (h==0){
      ofstream out;
      out.open(sampling+name1, ios::app);

      err_u=Error(glob_av[iu],glob_av2[iu],nblk);
      err_c=Error(glob_av[ic],glob_av2[ic],nblk);
      err_x=Error(glob_av[ix],glob_av2[ix],nblk);

      out << temp << "\t";
      out << glob_av[iu]/(double)nblk << "\t" << err_u << "\t";
      out << glob_av[ic]/(double)nblk << "\t" << err_c << "\t";
      out << glob_av[ix]/(double)nblk << "\t" << err_x << endl;

      out.close();
    }else{
      ofstream out;
      out.open(sampling+name2, ios::app);
      err_m=Error(glob_av[im],glob_av2[im],nblk);
      out << temp << "\t";
      out << glob_av[im]/(double)nblk << "\t" << err_m << endl;

    }
  }

}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
