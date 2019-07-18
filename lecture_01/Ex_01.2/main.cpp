#include "header.h"

using namespace std;

int main (int argc, char *argv[]){

	Input();

	ifstream numbers("numbers.in");									//input file containing 4 values of N, which defines the average variable S_N of the CLT

	string nomefile="histo_N=";
	string ext=".txt";

	double counter_std, counter_exp, counter_lrz;
	double std[10000];						//vectors containing 10^4 realizations of S_N
	double exp[10000];
	double lrz[10000];

	while ( !numbers.eof() ){
		ResetHisto();
		numbers >> N;
		FILE *histo;
		histo=fopen((nomefile + to_string(N) + ext).c_str(), "w");

		for (int i=0; i<10000; i++) {							//10^4 realizations of S_N
			counter_std=0;
			counter_exp=0;
			counter_lrz=0;
			for (int j=0; j<N; j++) {								//sum of x_j with j=1,...,N
				counter_std += rnd.Rannyu();
				counter_exp += rnd.Exponential(decay_rate);
				counter_lrz += rnd.Lorentz(width, median);
			}
			std[i]=counter_std/(double)N;
			exp[i]=counter_exp/(double)N;
			lrz[i]=counter_lrz/(double)N;
		}

		for (int i=0; i<10000; i++) {
			 fill_histo(nbin, 0., 1., h_std, std[i] );
			 fill_histo(nbin, 0., 6., h_exp, exp[i] );
			 fill_histo(nbin, -50., 50., h_lrz, lrz[i] );
		}

		double lbin_std=1./(double)nbin;
		double lbin_exp=6./(double)nbin;
		double lbin_lrz=100./(double)nbin;

		for (int i=0; i<nbin; i++) {								//write to "histo_N=N.txt"
		 		double x_std=lbin_std*(i+.5);
				double x_exp=lbin_exp*(i+.5);
				double x_lrz=-50+lbin_lrz*(i+.5);
			 fprintf(histo, "%d \t %f \t %d \t %f \t  %d \t %f \t %d \n", i, x_std, h_std[i], x_exp, h_exp[i], x_lrz, h_lrz[i]);
			 // col(0)=x_bin; col(1)=std_histogram; col(2)=exp_histogram; col(3)=lrz_histogram
		}

	}

   return 0;
}

void Input(void){

	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		 Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
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

}

void ResetHisto(void){
	for (int i=0; i<nbin; i++) {								//histograms inizialization
		h_std[i]=0;
		h_exp[i]=0;
		h_lrz[i]=0;
	}
}

void fill_histo (int nbin, double min, double max, int* h, double x){						// Attention: the vector h must be initialized to 0
	if (x<min || x>max){
		cout << x << " is not included in the range [" << min << "," << max <<"]" << endl;
		return;
	}
	double lbin=(max-min)/(double)nbin;
	for (int i=0; i<nbin; i++){
		if (x<min+(i+1)*lbin){ h[i]++; break; }
	}
	return;
}
