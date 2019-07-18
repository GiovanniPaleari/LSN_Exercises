#include "Map.h"
#include "random.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

Map::Map(){
	_Nc=0;
}

Map::Map(int N, int SorC, Random* rand, double l){
	_Nc=N;
  Pos r;
  if(SorC==0){      //cities are randomly placed inside a square
    for (int i=0; i<_Nc; i++){
      r.x=rand->Rannyu(-l/2.,l/2.);
      r.y=rand->Rannyu(-l/2.,l/2.);
      _cities.push_back(r);
    }
  }else{            //cities are randomly placed on a circumference
    for (int i=0; i<_Nc; i++){
      double theta = rand->Rannyu(0.,2*M_PI);
      r.x=l*cos(theta);
      r.y=l*sin(theta);
      _cities.push_back(r);
    }
  }

}

Map::~Map(){
	//delete[] _chromo;
}

double Map::GetX(int i){ return _cities.at(i).x; }
double Map::GetY(int i){ return _cities.at(i).y; }

double Map::sq_distance(int a,int b){

	double dist = pow((_cities.at(a).x-_cities.at(b).x),2.)+pow((_cities.at(a).y-_cities.at(b).y),2.);

	return dist;
}
