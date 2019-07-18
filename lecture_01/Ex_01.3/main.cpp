#include "header.h"

using namespace std;

int main (int argc, char *argv[]){

	FILE *pl;

	Input();

	int N_hit=0, N_thr=0;
	double x1, y1, alpha;
	double y2;
	double pi[N], pi2[N], cum_pi[N], cum_pi2[N], cum_error[N];

	for (int j=0; j<N; j++){		//reset averages
		pi[j]=0;
		pi2[j]=0;
		cum_pi[j]=0;
		cum_pi2[j]=0;
	}

	for (int j=0; j<N; j++){
		N_thr=0;
		N_hit=0;
		for(int i=0; i<block_size; i++){
			 //x1=rnd.Rannyu(0., size);
			 y1=rnd.Rannyu(0., size);		//throw of a needle, which is determined by x,y and direction
			 alpha=2*theta();

			 y2=y1+length*sin(alpha);
			 if(y2>y1 && ceil(y1)==floor(y2)){N_hit++;}
			 if(y1>y2 && ceil(y2)==floor(y1)){N_hit++;}

			 N_thr++;
		}
		pi[j]=2*length*N_thr/(double)(N_hit*dist);
		pi2[j]=pow(pi[j],2.);
	 }

	for (int j=0; j<N; j++){					//update cumulative averages
		 for(int i=0; i<j+1; i++){
			cum_pi[j] += pi[i];
			cum_pi2[j] += pi2[i];
		}
		cum_pi[j] /= (double)(j+1);
		cum_pi2[j] = cum_pi2[j]/(double)(j+1);
		cum_error[j] = error(cum_pi, cum_pi2, j);
	}

	pl=fopen("Plot_1.txt", "w");

	for (int j=0; j<N; j++){
		fprintf(pl, "%d \t %f \t %f \n", j, cum_pi[j], cum_error[j] );
	}
	fclose(pl);

   return 0;
}

void Input(void){

	ifstream ReadInput("Input.txt");
	ReadInput >> M >> N >> size >> dist >> length;
	block_size = M/N;
	cout << "Buffon experiment is performed with " << M <<" throws of the needle" << endl;
	cout << "The plane of the experimet is " << size << " long." << endl;
	cout << "The distance between lines is " << dist << endl;
	cout << "The needles are " << length << " long." << endl;

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

double theta (){				// random extraction of theta in [0,pi]
	double x=2, y=2;
	while (x*x+y*y>1) {
		x=rnd.Rannyu(-1.,1.);
		y=rnd.Rannyu();
	}

	return acos(x/sqrt(x*x+y*y));

}

double error(double* av, double* av2, int n){            // function for statistical uncertainty evaluation
	if (n==0)
		return 0;
	else
		return sqrt((av2[n]-pow(av[n],2.))/(double)n);
}
