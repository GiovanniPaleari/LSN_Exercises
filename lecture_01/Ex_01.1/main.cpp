#include "header.h"

using namespace std;

int main (int argc, char *argv[]){

	FILE *pl;

	Input();

	double ave[N], av2[N];                             // vectors of block (squared) averages
  double cum_ave[N], cum_av2[N], cum_error[N];       // vectors of cumulative (squared) averages and errors
	double counter=0;

	// f(x) = x

	for (int j=0; j<N; j++){
	   for(int i=0; i<block_size; i++){
      		counter += rnd.Rannyu();
		}
		ave[j]=counter/(double)block_size;	//update averages
		av2[j]=pow(ave[j],2.);
		counter=0;
   }

	for (int j=0; j<N; j++){
		cum_ave[j]=0;
		cum_av2[j]=0;
	   for(int i=0; i<j+1; i++){
			cum_ave[j] += ave[i];			//update cumulative averages
			cum_av2[j] += av2[i];
		}
		cum_ave[j] /= (double)(j+1);
		cum_av2[j] = cum_av2[j]/(double)(j+1);
		cum_error[j] = error(cum_ave, cum_av2, j);
	}

	pl=fopen("Plot_1.txt", "w");

	double c;
	for (int j=0; j<N; j++){
		c = cum_ave[j]-0.5;
		fprintf(pl, "%d \t %f \t %f \n", j, c, cum_error[j] );
	}
	fclose(pl);

	// f(x) = (x-1/2)^2

	double r;

	for (int j=0; j<N; j++){
		 for(int i=0; i<block_size; i++){
			 		r=rnd.Rannyu();
					counter += pow(r-0.5, 2.);
		}
		ave[j]=counter/(double)block_size;
		av2[j]=pow(ave[j],2.);
		counter=0;
	 }

	for (int j=0; j<N; j++){
		cum_ave[j]=0;
		cum_av2[j]=0;
		 for(int i=0; i<j+1; i++){
			cum_ave[j] += ave[i];
			cum_av2[j] += av2[i];
		}
		cum_ave[j] /= (double)(j+1);
		cum_av2[j] = cum_av2[j]/(double)(j+1);
		cum_error[j] = error(cum_ave, cum_av2, j);
	}

	pl=fopen("Plot_2.txt", "w");

	for (int j=0; j<N; j++){
		c = cum_ave[j] - 1/12.;
		fprintf(pl, "%d \t %f \t %f \n", j, c, cum_error[j] );
		}

  fclose(pl);

	int intervals[100];
	double chi_squared[100];
	double rr;
	counter=0;

	for(int k=0; k<100; k++){ intervals[k]=0; }								//initialize vector

	for(int k=0; k<100; k++){
		for (int i = 0; i < 10000; i++) {
			rr=rnd.Rannyu();
			intervals[int(rr*100)]++;		//fill histogram
		}
		counter=0;
		for (int l = 0; l < 100; l++) {
			counter += pow(intervals[l]-100., 2.)/100.;		//evaluate chi_squared
			intervals[l]=0;
		}
		chi_squared[k]=counter;
	}

	pl=fopen("Plot_Chi.txt", "w");

	for (int j=0; j<100; j++){
		fprintf(pl, "%d \t %f \n", j, chi_squared[j] );
		}

  fclose(pl);

  rnd.SaveSeed();
  return 0;
}

void Input(void){

	ifstream ReadInput("input.txt");
	ReadInput >> M;
	ReadInput >> N;

	block_size =M/N;

	cout << "The program estimates: " << endl;
	cout << "\t" << " - integral of x in (0,1)" << endl;
	cout << "\t" << " - integral of (x-0.5)^2 in (0,1)" << endl;
	cout << "using " << M << " throws divided in " << N << " blocks." << endl;

	ReadInput.close();


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

double error(double* av, double* av2, int n){            // function for statistical uncertainty evaluation
	if (n==0)
		return 0;
	else
		return sqrt((av2[n]-pow(av[n],2.))/(double)n);
}
