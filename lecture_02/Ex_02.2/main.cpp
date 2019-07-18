#include "header.h"

int main (int argc, char *argv[]){

	FILE *pl;

	Input();

	double xnow[M], ynow[M], znow[M];									// coordinates of current position for every simulation
	double current_sqdistance[M];											// squared distances of current random walkers in the M different simulations

//discrete 3D RW
	null(xnow,M);																			// simulations start at the origin
	null(ynow,M);
	null(znow,M);

	pl=fopen("Plot_lattice.txt", "w");

	for (int i=0; i<N; i++){
		for (int j=0; j<M; j++){
			double r=rnd.Rannyu(0.,3.);
			if(r<1.) xnow[j]=lattice_step(xnow[j]);
			else if(r<2. && r>1.) ynow[j]=lattice_step(ynow[j]);
			else if(r<3. && r>2.) znow[j]=lattice_step(znow[j]);
			current_sqdistance[j]=sq_distance(step_length*xnow[j], step_length*ynow[j], step_length*znow[j]);
		}
		fprintf(pl, "%d \t %f \t %f \n", i, sqrt(mean(current_sqdistance,M)), stddev(current_sqdistance,M));		//ATTENZIONE DA CAPIRE /(2*sqrt(mean(current_sqdistance,M)))
	}

	fclose(pl);

//continuum 3D RW
	null(xnow,M);																			// simulations start at the origin
	null(ynow,M);
	null(znow,M);

	double step_length=1.;
	pl=fopen("Plot_uniform.txt", "w");

	for (int i=0; i<N; i++){
		for (int j=0; j<M; j++){
			double phi=rnd.Rannyu(0.,2*M_PI);
			double theta=acos(1-2*rnd.Rannyu());			// p(theta)=.5*sin(theta)

			xnow[j]+=step_length*sin(theta)*cos(phi);
			ynow[j]+=step_length*sin(theta)*sin(phi);
			znow[j]+=step_length*cos(theta);
			current_sqdistance[j]=sq_distance(xnow[j], ynow[j], znow[j]);
		}
		fprintf(pl, "%d \t %f \t %f \n", i, sqrt(mean(current_sqdistance,M)), stddev(current_sqdistance,M));		//ATTENZIONE DA CAPIRE /(2*sqrt(mean(current_sqdistance,M)))
	}

	fclose(pl);

   rnd.SaveSeed();
   return 0;
}

void Input(){

	ifstream ReadInput("input.txt");
  ReadInput >> M;
  ReadInput >> N;

  cout << M << " simulations of a 3D random walk of " << N << " step on a cubic lattice and in the continuum." << endl;

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
	return;
}

double sq_distance(double x, double y, double z){ return (pow(x,2.)+pow(y,2.)+pow(z,2.));}

double mean(double* v, int N){
	double counter=0;
	for(int i=0; i<N; i++){counter+=v[i];}
	return counter/(double)N;
}

double stddev(double* v, int N){
	double counter=0;
	for(int i=0; i<N; i++){counter+=pow(v[i],2.);}
	return sqrt((counter/(double)N-pow(mean(v,N),2.))/(double)(N-1));
}

int lattice_step(int now){
	double r=rnd.Rannyu();
	if(r<.5){return now+step_length;}
	else{return now-step_length;}
	}

void null(double* v, int dim){
	for(int i=0; i<dim; i++)
		v[i]=0;
	return;
}
