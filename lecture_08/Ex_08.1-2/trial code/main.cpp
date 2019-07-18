#include "header.h"
#include "random.h"

int main (int argc, char *argv[]){

	Input();
	fiftypercent();
	Reset();

	for (int i=0; i<nblocks; i++){
		for (int j=0; j<block_size; j++){
			Metropolis();
			Measure();
		}
		Block(i);
	}

   return 0;
}

void Input(void){

	cout << "Single Quantum Particle in 1D confined by the external potential:        " << endl;
	cout << "V(x) = x^4 - 2.5 * x^2       " << endl;
	cout << "Variational Monte Carlo simulation             " << endl << endl;
	cout << "The trial wave function is the sum of 2 Gaussian functions:" << endl;
	cout << " exp(-(x-mu)^2/2*sigma^2) + exp(-(x+mu)^2/2*sigma^2) " << endl << endl;

//Random Inizialization
	int seed[4];
	int p1, p2;
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

//Read Input.dat
	input.open("input.dat");
	if (input.is_open()){
		 input >> step >> sigma >> mu >> nblocks >> nstep >> xnow;
	} else cerr << "PROBLEM: Unable to open Input.dat" << endl;
	input.close();

	cout << "Number of Monte Carlo steps: " << nstep << endl;
	cout <<"Number of blocks: " << nblocks << endl << endl;
	cout << "Initial position: x = " << xnow << endl << endl;

	if (nstep%nblocks != 0) cout << "*** ACHTUNG: nstep non è multiplo di nblocks ***" << endl;
	block_size = nstep/nblocks;

	ave=0;
	cum_ave=0;		//Reset averages
	cum_av2=0;
	istep=0;

	return;
}

void Reset(void){

	ave=0;
	cum_ave=0;		//Reset averages
	cum_av2=0;
	istep=0;

	return;
}

void Metropolis(void){

		xtry = xnow + rnd.Rannyu(-step, step);

	acc = min(1.,pow(FTrial(xtry)/FTrial(xnow),2.));

	double rand = rnd.Rannyu();
	if(rand<=acc){
		xnow=xtry;
	}

	istep++;

	return;
}

void Measure(void){
	double Vpsi = (pow(xnow,4.) - 2.5*pow(xnow,2.))*FTrial(xnow);
  double Tpsi = -1/pow(sigma,2.)*FTrial(xnow) + pow(xnow-mu,2.)/pow(sigma,4.)*exp(-.5*pow((xnow-mu)/sigma,2.)) + pow(xnow+mu,2.)/pow(sigma,4.)*exp(-.5*pow((xnow+mu)/sigma,2.));
  Hpsi = (Tpsi + Vpsi)/FTrial(xnow);

	ave += Hpsi;

	/*if (istep%iprint==0){
		if(istep==100) output.open("actual_pos_dist.txt");
		else output.open("actual_pos_dist.txt", ios::app);
		output << istep << "\t" << xnow << "\t" << ynow << "\t" << znow << "\t" << rnow << endl;
		output.close();
	}*/
	return;
}

void Block(int iblock){			// Data blocking
	ave /= (double)block_size;
	cum_ave += ave;
	cum_av2 += ave*ave;

	if(iblock==0) output.open("energy.txt");
	else output.open("energy.txt", ios::app);

	double err = error(cum_ave/(double)(iblock+1), cum_av2/(double)(iblock+1), iblock);

	output << iblock << "\t" << cum_ave/(double)(iblock+1) << "\t" << err << endl;
	output.close();

	ave=0;
	return;
}

double FTrial(double x){
  double f = exp(-pow(x-mu,2.)/(2*pow(sigma,2.)))+exp(-pow(x+mu,2.)/(2*pow(sigma,2.)));
  return f;
}

void fiftypercent(void){

		double q=0;
		while (q>.55 || q<.45){
			int counter=0;
			for(int k=0; k<1000; k++){
				Metropolis();
				if(xnow==xtry) counter++;		//contatore delle mosse accettate
			}
			q=counter/1000.;
			step += .1;
		}
cout << "In order tu fulfil 50% rule step is fixed to: " << step << endl;


	return;
}

double error(double av, double av2, int n){            // function for statistical uncertainty evaluation
	if (n==0)
		return 0;
	else
		return sqrt((av2-pow(av,2.))/(double)n);
}
