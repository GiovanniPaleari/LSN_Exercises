#include "header.h"

using namespace std;

//double error (double*, double*, int);

int main (int argc, char *argv[]){

	FILE *pl;

	double ave[N], av2[N];                             // vectors of block (squared) averages
  double cum_ave[N], cum_av2[N], cum_error[N];       // vectors of cumulative (squared) averages and errors
	double counter=0;

	// f(x) = pi/2*cos(pi*x/2)

	//Uniform distribution

	for (int j=0; j<N; j++){
	   for(int i=0; i<block_size; i++){
      		counter += M_PI/2.*cos(M_PI*rnd.Rannyu()/2.);
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

	pl=fopen("Plot_uniform.txt", "w");

	double c;
	for (int j=0; j<N; j++){
		c = cum_ave[j]-1;
		fprintf(pl, "%d \t %f \t %f \n", j, c, cum_error[j] );
	}
	fclose(pl);

	//Importance Sampling -- p(x)=2*(1-x)

	counter=0;
	for (int j=0; j<N; j++){
	   for(int i=0; i<block_size; i++){
			 double x=1-sqrt(1-rnd.Rannyu());				 						//sampling p(x) through the inverse of the cumulative
			 counter += M_PI/4.*cos(M_PI*x/2.)/(double)(1-x);		//evaluate f(x)/p(x) in points sampled from p(x)
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

	pl=fopen("Plot_imp_samp.txt", "w");


	for (int j=0; j<N; j++){
		c = cum_ave[j]-1;
		fprintf(pl, "%d \t %f \t %f \n", j, c, cum_error[j] );
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

  cout << "The program estimates the integral of pi/2*cos(pi*x/2) in (0,1)" << endl;
  cout << "using " << M << " throws divided in " << N << " blocks." << endl;
	cout << "The calculation is performed in 2 ways:" << endl;
	cout << "\t - by sampling a uniform distribution in [0,1]" << endl;
	cout << "\t - using importance sampling" << endl;

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
