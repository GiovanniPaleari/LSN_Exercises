#include "header.h"
#include "random.h"

int main (int argc, char *argv[]){

	Input();
	Reset();


	for(int iblock=0; iblock<N; iblock++){
		call(iblock,block_size,1,"Call_1step.txt");
		put(iblock,block_size,1,"Put_1step.txt");
	}

	Reset();

	for(int iblock=0; iblock<N; iblock++){
		call(iblock,block_size,100,"Call_100step.txt");
		put(iblock,block_size,100,"Put_100step.txt");
	}

   return 0;
}

void Reset(void){		//reset cumulative averages
	cum_ave[iput]=0;
	cum_av2[iput]=0;
	cum_ave[icall]=0;
	cum_av2[icall]=0;
}


void Input(void){

	input.open("Input.dat");
	if (input.is_open()){
		 input >> M >> N >> S_0 >> T >> K >> r >> sigma;
	} else cerr << "PROBLEM: Unable to open Input.dat" << endl;
	input.close();
	block_size = M/N;

	cout << "Plain Vanilla Option pricing" << endl;
	cout << "The simulation uses the following parameters:" << endl;
	cout << "\t - asset price at t=0: S_0 = " << S_0 << endl;
	cout << "\t - delivery time: T = " << T << endl;
	cout << "\t - strike price: K = " << K << endl;
	cout << "\t - risk-free interest rate: r = " << r << endl;
	cout << "\t - volatility: sigma = " << sigma << endl << endl;
	cout << "The program performs " << M << " simulations." << endl;

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
	return;
}

double error(double av, double av2, int n){            // function for statistical uncertainty evaluation
	if (n==0)
		return 0;
	else
		return sqrt((av2-pow(av,2.))/(double)n);
}


void call(int iblock, int block_size, int nstep, string nomefile){
	ave[icall]=0;
	double S_t;
	double delta_t=T/(double)nstep;
	for (int i=0; i<block_size; i++){
		S_t=S_0;
		for (int j=0; j<nstep; j++){
			S_t=S_t*exp((r-.5*pow(sigma,2.))*delta_t + sigma*rnd.Gauss(0.,1.)*sqrt(delta_t));
		}
		ave[icall] += exp(-r*T)*fmax(0.,(S_t-K));
	}
	ave[icall]=ave[icall]/(double)block_size;
	cum_ave[icall] += ave[icall];
	cum_av2[icall] += pow(ave[icall], 2.);

	if (iblock==0){
		out=fopen(nomefile.c_str(), "w");										//open file in write mode
		fprintf(out, "%d \t %f \t %f \n", iblock, cum_ave[icall]/(double)(iblock+1), error(cum_ave[icall]/(double)(iblock+1), cum_av2[icall]/(double)(iblock+1), iblock));
	}else{
		out=fopen(nomefile.c_str(), "a");										//open file in append mode
		fprintf(out, "%d \t %f \t %f \n", iblock, cum_ave[icall]/(double)(iblock+1), error(cum_ave[icall]/(double)(iblock+1), cum_av2[icall]/(double)(iblock+1), iblock));
	}
	fclose(out);
	return;
}

void put(int iblock, int block_size, int nstep, string nomefile){
	ave[iput]=0;
	double S_t;
	double delta_t=T/(double)nstep;
	for (int i=0; i<block_size; i++){
		S_t=S_0;
		for (int j=0; j<nstep; j++){
			S_t=S_t*exp((r-.5*pow(sigma,2.))*delta_t + sigma*rnd.Gauss(0.,1.)*sqrt(delta_t));
		}
		ave[iput] += exp(-r*T)*fmax(0.,(K-S_t));
	}
	ave[iput]=ave[iput]/(double)block_size;
	cum_ave[iput] += ave[iput];
	cum_av2[iput] += pow(ave[iput], 2.);

	if (iblock==0){
		out=fopen(nomefile.c_str(), "w");										//open file in write mode
		fprintf(out, "%d \t %f \t %f \n", iblock, cum_ave[iput]/(double)(iblock+1), error(cum_ave[iput]/(double)(iblock+1), cum_av2[iput]/(double)(iblock+1), iblock));
	}else{
		out=fopen(nomefile.c_str(), "a");										//open file in append mode
		fprintf(out, "%d \t %f \t %f \n", iblock, cum_ave[iput]/(double)(iblock+1), error(cum_ave[iput]/(double)(iblock+1), cum_av2[iput]/(double)(iblock+1), iblock));
	}
	fclose(out);
	return;
}
