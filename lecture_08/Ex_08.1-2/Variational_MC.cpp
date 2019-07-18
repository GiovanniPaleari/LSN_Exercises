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
#include <cmath>
#include <iomanip>
#include "Variational_MC.h"

using namespace std;

int main()
{
  Input(); //Inizialization
  fiftypercent();

/*
  for (int i=0; i<1000; i++){   //Equilibration
    Move();
  }*/

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }

  PrintHisto();

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Single Quantum Particle in 1D confined by the external potential:        " << endl;
  cout << "V(x) = x^4 - 2.5 * x^2       " << endl;
  cout << "Variational Monte Carlo simulation             " << endl << endl;
  cout << "The trial wave function is the sum of 2 Gaussian functions:" << endl;
  cout << " exp(-(x-mu)^2/2*sigma^2) + exp(-(x+mu)^2/2*sigma^2) " << endl << endl;

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

  ReadInput >> jump;

  ReadInput >> sigma;
  cout << "Sigma = " << sigma << endl;

  ReadInput >> mu;
  cout << "Mu = " << mu << endl;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> xin;
  xnow = xin;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << jump << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  cout << "Initial Position: x = " << xin << endl << endl;
  ReadInput.close();

  //Prepare histogram from -2.5 to 2.5
  bin_size=5./(double)nbin;

  for (int i=0; i<100; i++) histo[i]=0;


}


void Move(void)
{

  double r = rnd.Rannyu(-1*jump,jump);

  xtry = xnow + r;

  double rr = rnd.Rannyu();
  double rapp = pow((FTrial(xtry)/FTrial(xnow)),2.);

  if (rr<rapp){
    xnow=xtry;
    accepted = accepted + 1.;
  }

  //fill histogram
  double bin = (xnow+2.5)/5.*100;   //rescale bin from [-2.5,2.5] to [0,5] and then to [0,100]

  if (xnow>-2.5 && xnow<2.5) histo[(int)floor(bin)]++;

  attempted = attempted + 1.;

}

void Measure()
{
  double Vpsi = (pow(xnow,4.) - 2.5*pow(xnow,2.))*FTrial(xnow);
  double Tpsi = -0.5*(-1/pow(sigma,2.)*FTrial(xnow) + pow(xnow-mu,2.)/pow(sigma,4.)*exp(-.5*pow((xnow-mu)/sigma,2.)) + pow(xnow+mu,2.)/pow(sigma,4.)*exp(-.5*pow((xnow+mu)/sigma,2.)));
  Hpsi = (Tpsi + Vpsi)/FTrial(xnow);

}


void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
     glob_av = 0;
     glob_av2 = 0;
   }


   blk_av = 0;
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{
  blk_av = blk_av + Hpsi;
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

  ofstream energy;
  string s = to_string(sigma);
  string m = to_string(mu);
  s.erase (3, std::string::npos);
  m.erase (3, std::string::npos);
  //const int wd=12;


  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;

  if (iblk==1) energy.open("Energy_sigma"+s+"_mu"+m+".txt");
  else energy.open("Energy_sigma"+s+"_mu"+m+".txt",ios::app);

  stima = blk_av/blk_norm;
  glob_av += stima;
  glob_av2 += stima*stima;
  err=Error(glob_av,glob_av2,iblk);

  energy << iblk <<  "\t" << glob_av/(double)iblk <<  "\t" << err << endl;


  cout << "----------------------------" << endl << endl;

  energy.close();
}

double Error(double sum, double sum2, int iblk)
{
  if (iblk==1)
    return 0;
  else
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));

}

double FTrial(double x){
  double f = exp(-pow(x-mu,2.)/(2*pow(sigma,2.)))+exp(-pow(x+mu,2.)/(2*pow(sigma,2.)));
  return f;
}

void fiftypercent(void){
    jump = 0.1;
		double q=0;
		while (q>.6 || q<.4){
      jump += 0.1;
      attempted = 0;
      accepted = 0;
			for(int k=0; k<10000; k++){
				Move();
			}
			q = accepted/attempted;

      //cout <<  jump << " " << q << endl;

    }

cout << "In order tu fulfil 50% rule step is fixed to: " << jump << endl;
}

void PrintHisto(void){

  ofstream h;
  h.open("histo.txt");

  int counter=0;
  for(int i=0; i<nbin; i++){
      counter += histo[i];
  }

  for(int i=0; i<nbin; i++){
    double xbin = i/100.*5.-2.5;
    h << xbin << "\t" << histo[i]/(counter*bin_size) << endl;
  }

  h.close();

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
