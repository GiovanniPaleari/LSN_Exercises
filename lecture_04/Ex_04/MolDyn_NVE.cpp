/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

int main(){
  Input();                //Inizialization
  int nconf = 1;

  for(int istep=1; istep <= nstep; ++istep){
     Move();              //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();        //Properties measurement
  //    ConfXYZ(nconf);   //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf += 1;
     }
  }
  ConfFinal();         //Write final configuration to restart

  Deallocate();

  return 0;
}

void rand_inizialization(void){
	int seed[4];
	int p1, p2;
  ifstream input;
	input.open("Primes");
	if (input.is_open()){
		 input >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	input.close();

	input.open("seed.in");
	string property;
	if (input.is_open()){
		 while ( !input.eof() ){
				input >> property;
				if( property == "RANDOMSEED" ){
					 input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
					 rnd.SetRandom(seed,p1,p2);
				}
		 }
		 input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	 rnd.SaveSeed();
	return;
}


void Input(void){       //Prepare all the stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  x = new double [npart];
  y = new double [npart];
  z = new double [npart];
  xold = new double [npart];
  yold = new double [npart];
  zold = new double [npart];
  vx = new double [npart];
  vy = new double [npart];
  vz = new double [npart];

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;

  ReadInput >> restart;
  if(restart==true)  cout << "Re-start option ENABLED" << endl << endl;
  else  cout << "Re-start option DISABLED" << endl << endl;

  ReadInput >> nblocks;
  if((nstep/10)%nblocks != 0) cout << "*** WARNING: the number of blocks is not a multiple of the number of steps." << endl << endl;
  //1 measurement is performed every 10 MD steps

  ReadInput.close();

//Reset (cumulative) averages
  ave_epot=0;
  ave_ekin=0;
  ave_pr=0;
  cum_av2_epot=0;
  cum_av2_ekin=0;
  cum_av2_pr=0;
  cum_ave_epot=0;
  cum_ave_ekin=0;
  cum_ave_pr=0;

//Reset the counter of the measurements performed on the system
  imeasure=0;

//block size inizialization. Remember: nstep MD steps --> nstep/10 MD measurements
  block_size=nstep/(10.*(double)nblocks);

  rand_inizialization();

  if(restart==true){

    //Read initial ACTUAL configuration
      cout << "Read initial actual configuration from file old.0 " << endl << endl;
      ReadConf.open("old.0");
      for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
      }
      ReadConf.close();

    //Read initial OLD configuration
    cout << "Read initial old configuration from file old.final " << endl << endl;
    ifstream ReadOldConf;
    ReadOldConf.open("old.final");
    for (int i=0; i<npart; ++i){
      ReadOldConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadOldConf.close();

    CorrectOldConf();

  }else if(restart==false){

    //Read initial configuration
      cout << "Read initial configuration from file config.0 " << endl << endl;
      ReadConf.open("config.0");
      for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
      }
      ReadConf.close();

//Prepare initial velocities if restart option is disabled
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = rnd.Rannyu() - 0.5;
      vy[i] = rnd.Rannyu() - 0.5;
      vz[i] = rnd.Rannyu() - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }

   return;
}

void CorrectOldConf(void){
  double xnew[npart], ynew[npart], znew[npart], fx[npart], fy[npart], fz[npart];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Compute r(t+dt) with 1 step of Verlet algorithm and compute v(t+dt/2)

    xnew[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew[i] - x[i])/( delta);
    vy[i] = Pbc(ynew[i] - y[i])/( delta);
    vz[i] = Pbc(znew[i] - z[i])/( delta);
  }

  //Compute the effective temperature
  double sumv2=0.;
  for (int i=0; i<npart; ++i)  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];

  sumv2 /= (double)npart;

  //Rescale the velocities in order to match the desired temperature
  double eff_temp= sumv2/3.;
  double fs=sqrt(temp/eff_temp);

  cout << "Target Temperature = " << temp << endl;
  cout << "Actual Temperature = " << eff_temp << endl;
  cout << "Scale Factor = " << fs << endl;

  for(int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;
  }

  //use v rescaled to estimate a novel old spatial configuration

  for(int i=0; i<npart; ++i){
    xold[i]=Pbc(xnew[i]-delta*vx[i]);
    yold[i]=Pbc(ynew[i]-delta*vy[i]);
    zold[i]=Pbc(znew[i]-delta*vz[i]);
    x[i]=xnew[i];
    y[i]=ynew[i];
    z[i]=znew[i];

  }

  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[npart], fy[npart], fz[npart];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(){ //Properties measurement

  imeasure++; //this is the (imeasure)th measurement performed on the system
  int iblock=imeasure/block_size;

  //int bin;
  double v, t, pr, vij, prij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pr;
  ofstream Ave_Epot, Ave_Ekin, Ave_Etot, Ave_Temp, Ave_Pr;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pr.open("output_pr.dat",ios::app);

  Ave_Epot.open("ave_epot.out",ios::app);
  Ave_Ekin.open("ave_ekin.out",ios::app);
  Ave_Temp.open("ave_temp.out",ios::app);
  Ave_Etot.open("ave_etot.out",ios::app);
  Ave_Pr.open("ave_pr.out",ios::app);

  if((imeasure-1)%block_size==0){   //this condition becomes true when a new block is starting
    ave_epot=0; //reset averages
    ave_ekin=0;
    ave_pr=0;
  }


  v = 0.0; //reset observables
  t = 0.0;
  pr= 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       prij = 48.*(1./pow(dr,12) - .5/pow(dr,6));

//Potential energy
       v += vij;
//Pressure
       pr += prij;
     }
    }
  }
  ave_epot += v;

//Kinetic energy
  for (int i=0; i<npart; ++i){
    t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }
  ave_ekin += t;



    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy
    stima_pr = rho*stima_temp + pr*rho/(double)(3*npart); //Pressure

    ave_pr += stima_pr;

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pr << stima_pr << endl;

    if(imeasure%block_size==0){   //this condition is true if this is the last measurement of the current block
      //Potential energy
      ave_epot /= (double)block_size*npart;
      cum_ave_epot += ave_epot;
      cum_av2_epot += pow(ave_epot,2.);

      double er_pot = error(cum_ave_epot/(double)(iblock), cum_av2_epot/(double)(iblock), iblock-1);
      Ave_Epot << cum_ave_epot/(double)(iblock) << "\t" << er_pot  << endl;

      //Kinetic energy
      ave_ekin /= (double)block_size*npart;
      cum_ave_ekin += ave_ekin;
      cum_av2_ekin += pow(ave_ekin,2.);

      double er_kin = error(cum_ave_ekin/(double)(iblock), cum_av2_ekin/(double)(iblock), iblock-1);
      Ave_Ekin << cum_ave_ekin/(double)(iblock) << "\t" << er_kin << endl;

      //Temperature
      Ave_Temp << (2.0 / 3.0) * cum_ave_ekin/(double)(iblock) << "\t" << (2.0 / 3.0) * er_kin << endl;

      //Total Energy
      Ave_Etot << cum_ave_ekin/(double)(iblock) + cum_ave_epot/(double)(iblock) << "\t" << sqrt(pow(er_kin,2.)+pow(er_pot,2.)) << endl;

      //Pressure
      ave_pr /= (double)block_size;
      cum_ave_pr += ave_pr;
      cum_av2_pr += pow(ave_pr,2.);

      double er_pr = error(cum_ave_pr/(double)(iblock), cum_av2_pr/(double)(iblock), iblock-1);
      Ave_Pr << cum_ave_pr/(double)(iblock) << "\t" << er_pr  << endl;
      }

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pr.close();

    Ave_Epot.close();
    Ave_Ekin.close();
    Ave_Temp.close();
    Ave_Etot.close();
    Ave_Pr.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

 //Ex 04.1.1: save the actual (old.0) and the old (old.final) spatial configurations in order to restart the simulation

    ofstream WriteOldConf;
    WriteConf.open("old.0");
    WriteOldConf.open("old.final");
    for (int i=0; i<npart; ++i){
      WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
      WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }
    WriteConf.close();
    WriteOldConf.close();


  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double error(double av, double av2, int n){            // function for statistical uncertainty evaluation
	if (n==0)
		return 0;
	else
		return sqrt((av2-pow(av,2.))/(double)n);
}

void Deallocate(void){
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] xold;
  delete[] yold;
  delete[] zold;
  delete[] vx;
  delete[] vy;
  delete[] vz;
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
