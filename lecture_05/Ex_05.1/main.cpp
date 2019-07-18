#include "header.h"
#include "random.h"

int main (int argc, char *argv[]){

	Input();
	//fiftypercent();
	//Equilibration(500);
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

	cout << "Sampling of wave functions of H-atom via Metropolis Algorithm" << endl << endl;
	cout << "The program uses Bohr radius units for distances" << endl << endl;

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
	input.open("Input.dat");
	if (input.is_open()){
		 input >> xnow >> ynow >> znow >> nstep >> nblocks >> T >> step >> n;
	} else cerr << "PROBLEM: Unable to open Input.dat" << endl;
	input.close();

	if (T==0) cout << "Uniform transition probability in (-" << step << "," << step << ")" << endl << endl;
	else if (T==1) cout << "Multivariate normal transition probability" << endl << endl;

	if (n==0) cout << "Sampling of 1s wave function" << endl << endl;
	else if (n==1) cout << "Sampling of 2p wave function" << endl << endl;

	cout << "Number of Monte Carlo steps: " << nstep << endl;
	cout <<"Number of blocks: " << nblocks << endl << endl;
	cout << "Initial position: (" << xnow << "," << ynow << "," << znow << ")" << endl << endl;

	if (nstep%nblocks != 0) cout << "*** ACHTUNG: nstep non è multiplo di nblocks ***" << endl;
	block_size = nstep/nblocks;

	ave=0;
	cum_ave=0;		//Reset averages
	cum_av2=0;
	istep=0;

	return;
}

void Equilibration(int eq){
	for (int i=0; i<eq; i++) Metropolis();
}

void Reset(void){

	ave=0;
	cum_ave=0;		//Reset averages
	cum_av2=0;
	istep=0;

	return;
}

void Metropolis(void){

	if (T==0){				//Uniform transition probability
		xtry = xnow + rnd.Rannyu(-step, step);
		ytry = ynow + rnd.Rannyu(-step, step);
		ztry = znow + rnd.Rannyu(-step, step);
	}else if(T==1){		//Multivariate normal transition probability
		xtry = xnow + rnd.Gauss(0.,step);
		ytry = ynow + rnd.Gauss(0.,step);
		ztry = znow + rnd.Gauss(0.,step);
	}

	if (n==0) acc = min(1.,pow(psi_1s(xtry, ytry, ztry)/psi_1s(xnow, ynow, znow),2.));				//1s wave function
	else if (n==1) acc = min(1.,pow(psi_2p(xtry, ytry, ztry)/psi_2p(xnow, ynow, znow),2.));		//2p wave function

	double rand = rnd.Rannyu();
	if(rand<=acc){
		xnow=xtry;
		ynow=ytry;
		znow=ztry;
	}

	istep++;

	return;
}

void Measure(void){
	rnow = sqrt(xnow*xnow + ynow*ynow + znow*znow);
	ave += rnow;

	if (istep%iprint==0){
		if(istep==100) output.open("actual_pos_dist.txt");
		else output.open("actual_pos_dist.txt", ios::app);
		output << istep << "\t" << xnow << "\t" << ynow << "\t" << znow << "\t" << rnow << endl;
		output.close();
	}
	return;
}

void Block(int iblock){			// Data blocking
	ave /= (double)block_size;
	cum_ave += ave;
	cum_av2 += ave*ave;

	if(iblock==0) output.open("mean_radius.txt");
	else output.open("mean_radius.txt", ios::app);

	double err = error(cum_ave/(double)(iblock+1), cum_av2/(double)(iblock+1), iblock);

	output << iblock << "\t" << cum_ave/(double)(iblock+1) << "\t" << err << endl;
	output.close();

	cout << "Block: " << iblock << "/" << nblocks-1 << endl;

	ave=0;
	return;
}

double psi_1s(double x, double y, double z){
	double r2=x*x+y*y+z*z;
	return exp(-sqrt(r2))/sqrt(M_PI);
}

double psi_2p(double x, double y, double z){
	double r2=x*x+y*y+z*z;
	return 1./8.*sqrt(2./M_PI)*sqrt(r2)*exp(-sqrt(r2)/2.)*z/sqrt(r2);
}

void fiftypercent(void){
	step=0.1;
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
