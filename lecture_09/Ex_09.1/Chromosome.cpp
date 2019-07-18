#include "Chromosome.h"
#include "random.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;

Chromosome::Chromosome(){
	_N=0;
}

Chromosome::Chromosome(int N, int mut, double pcross, double* pmut, Random* rnd){
	_N=N;
	_mut=mut;
	_rnd=rnd;

	_p_mutation = new double[mut];
	for(int i=0; i<mut; i++) _p_mutation[i]=pmut[i];

	_p_crossover = pcross;

	for (int i=0; i<_N; i++) _chromo.push_back(i);

	random_shuffle ( _chromo.begin(), _chromo.end() );
}

Chromosome::Chromosome(Chromosome* bro){
	_N=bro->GetN();
	_mut=bro->GetMut();
	_rnd=bro->GetRnd();

	_p_mutation = new double[_mut];


	for(int i=0; i<_mut; i++) _p_mutation[i]=bro->GetP(i);

	_p_crossover=bro->GetPcross();

	for (int i=0; i<_N; i++) _chromo.push_back(bro->GetValue(i));

}

Chromosome::~Chromosome(){
}

int Chromosome::GetValue(int pos){ return _chromo.at(pos); }

void Chromosome::SetValue(int pos, int value){ _chromo.at(pos)=value; }

void Chromosome::ChangeProb(){
	for (int i=0; i<_mut; i++) _p_mutation[i]=_rnd->Rannyu(0.,0.1);
  _p_crossover=_rnd->Rannyu(0.5,1.);
}

int Chromosome::GetN(){ return _N; }
int Chromosome::GetMut(){ return _mut; }
double Chromosome::GetP(int i){ return _p_mutation[i]; }
double Chromosome::GetPcross(){ return _p_crossover; }
Random* Chromosome::GetRnd(){ return _rnd; }
std::vector<int>::iterator Chromosome::begin(){ return _chromo.begin(); }

Chromosome Chromosome::operator=(const Chromosome& c){
	_N=c._N;
	_mut=c._mut;
	for(int i=0; i<_mut; i++) _p_mutation[i]=c._p_mutation[i];
	_p_crossover=	c._p_crossover;
	for (int i=0; i<_N; i++) _chromo.at(i)=c._chromo.at(i);
	_rnd=c._rnd;

	return *this;
}

void Chromosome::PrintAdjList(string nomefile){
	ofstream out;
	out.open(nomefile);
	out << _chromo.at(0) << "\t" << _chromo.at(_N-1) << endl;
	for(int i=0; i<_N-1; i++){
		out << _chromo.at(i) << "\t" << _chromo.at(i+1) << endl;
		//if(i%2!=0) out << endl;
	}
	out.close();
}

int Chromosome::CheckBonds(void){

	//the salesman must visit one and only one time every city

	vector <int> tmp(_N);

	copy(_chromo.begin(), _chromo.end(), tmp.begin());

	sort(tmp.begin(), tmp.end());

	if (distance(tmp.begin(),unique(tmp.begin(), tmp.end()))!=_N){
		cout << "Error: bond not fulfilled " << endl;
		return -1;
	}

	return 0;

}

double Chromosome::Fitness_L2(Map cities){

	double path=cities.sq_distance(_chromo.front(),_chromo.back());

	for(int i=0; i<_N-1; i++){
		path += cities.sq_distance(_chromo.at(i),_chromo.at(i+1));
	}

	return path;
}

double Chromosome::Fitness_L1(Map cities){

	double path=pow(cities.sq_distance(_chromo.front(),_chromo.back()),.5);

	for(int i=0; i<_N-1; i++)  path += pow(cities.sq_distance(_chromo.at(i),_chromo.at(i+1)),.5);

	return path;
}

void Chromosome::PairPermutation(){
	double rand1 = _rnd->Rannyu();
	//cout << rand << endl;
	if (rand1<_p_mutation[0]){
		double rand2 = _rnd->Rannyu(0.,_N);
		double rand3 = _rnd->Rannyu(0.,_N);
		swap(_chromo.at(floor(rand2)),_chromo.at(floor(rand3)));
	}
}

void Chromosome::Shift(){
	double rand1 = _rnd->Rannyu();
	if (rand1<_p_mutation[1]){
		double rand2 = _rnd->Rannyu(0.,_N);
		int n = floor(rand2);
		rotate(_chromo.begin(),_chromo.end()-n,_chromo.end());
	}
}

void Chromosome::RigidShift(){
	//rand1 is a random number in (0,1)
	//rand2 is a random number in (0,ncities-n-m)
	//rigid shift of +n positions of m contiguous cities

	double rand1 = _rnd->Rannyu();

	double rand2 = _rnd->Rannyu(0.,_N/(double)2.);		//m
	int m = floor(rand2);

	double rand3 = _rnd->Rannyu(0.,_N/(double)2.);				//n
	int n = floor(rand3);

	double rand4 = _rnd->Rannyu(0.,_N-n-m);

	int pos = floor(rand4);

	if (pos+m+n >= _N){ cout << "Error: Size exceeded" << endl; return; }

	if (rand1<_p_mutation[2]){
		vector<int> tmp(n);
		copy(_chromo.begin()+pos+m,_chromo.begin()+pos+m+n,tmp.begin());
		move_backward(_chromo.begin()+pos,_chromo.begin()+pos+m, _chromo.begin()+pos+m+n);
		copy(tmp.begin(),tmp.end(),_chromo.begin()+pos);
	}

}

void Chromosome::SwapRange(){

	double rand0 = _rnd->Rannyu();
	double rand1 = _rnd->Rannyu(0,_N/(double)2.);
	int m = floor(rand1);
	if (rand0<_p_mutation[3]){
		double rand2 = _rnd->Rannyu(0.,floor(_N/2.)-m);
		int pos2 = floor(rand2);
		double rand3 = _rnd->Rannyu(pos2+m,_N-m);
		int pos3 = floor(rand3);
		swap_ranges(_chromo.begin()+pos2, _chromo.begin()+pos2+m, _chromo.begin()+pos3);
	}

}

void Chromosome::Reverse(){

	double rand0 = _rnd->Rannyu();
	double rand1 = _rnd->Rannyu(0,_N);
	int m = floor(rand1);
	if (rand0<_p_mutation[4]){
		double rand2 = _rnd->Rannyu(0.,_N-m);
		int pos = floor(rand2);
		reverse(_chromo.begin()+pos,_chromo.begin()+pos+m);
	}

}
