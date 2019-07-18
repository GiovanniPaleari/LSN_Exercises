#ifndef __Chromosome_h__
#define __Chromosome_h__
//#include "GA_tsp.h"

#include "Map.h"
#include <vector>
#include <string>
#include <fstream>

#include "random.h"

class Chromosome{

public:
	Chromosome();
	Chromosome(int, int, double, double*, Random*);
	Chromosome(Chromosome*);
	~Chromosome();

	int CheckBonds(void);
	int GetValue(int);
	void SetValue(int, int);
	void ChangeProb();
	int GetN();
	int GetMut();
	double GetP(int);
	double GetPcross();
	Random* GetRnd();
	std::vector<int>::iterator begin();

	Chromosome operator=(const Chromosome&);

	void PrintAdjList(std::string);

	double Fitness_L2(Map);
	double Fitness_L1(Map);

	//mutation operators
	void PairPermutation();
	void Shift();
	void RigidShift();
	void SwapRange();
	void Reverse();


protected:
	int _N;							//number of genes
	int _mut;						//number of mutation operators
	std::vector<int> _chromo;
	double* _p_mutation;
	double _p_crossover;
	Random* _rnd;

};

#endif
