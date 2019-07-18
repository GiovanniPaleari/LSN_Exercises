#include "GA_tsp.h"

using namespace std;

int main()
{
  Input();
  Initialize();
  best.open("BestPath.txt");

  for (int i=0; i<generations; i++){
    cout << "----------------------------------" << endl;

    cout  << "Generation " << i << endl << endl;
    for (int j=0; j<npop; j=j+2){
      int u=Selection();
      int v;
      do{ v = Selection(); }while(v==u);    //select the second chromosome different from the first
      newgen.at(j)=population.at(u);
      newgen.at(j+1)=population.at(v);

      Crossover(population.at(u), population.at(v), &newgen.at(j), &newgen.at(j+1));

      //check
      newgen.at(j).CheckBonds();
      newgen.at(j+1).CheckBonds();

      //mutations
      newgen.at(j).PairPermutation();
      newgen.at(j+1).PairPermutation();
    	newgen.at(j).Shift();
      newgen.at(j+1).Shift();
    	newgen.at(j).RigidShift();
      newgen.at(j+1).RigidShift();
    	newgen.at(j).SwapRange();
      newgen.at(j+1).SwapRange();
    	newgen.at(j).Reverse();
      newgen.at(j+1).Reverse();

    }

    for(int t=0; t<npop; t++) population.at(t) = newgen.at(t);

    sort(population.begin(),population.end(),order);

    double counter=0;
    for(int t=0; t<npop/2; t++) counter += population.at(t).Fitness_L1(*cities);

    best << i << "\t" << counter*2/npop << "\t" << population.at(0).Fitness_L1(*cities) << endl;

    cout << "Best Path: " << population.at(0).Fitness_L1(*cities) << endl;
  }

  population.at(0).PrintAdjList("AdjList.txt");   //print adjacency list of the best path at the end of the algorithm

  best.close();
  Deallocate();
  return 0;
}

void Input(void)
{

  cout << "Genetic Algorithm to optimize the Traveling Salesman Problem" << endl << endl;
  cout << "Chromosomes are the possible closed path which visit each city only once" << endl << endl;

  //Read seed for random numbers
  rnd = new Random;

  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd->SetRandom(seed,p1,p2);
  input.close();

  ifstream ReadInput;


  ReadInput.open("input.dat");

  ReadInput >> npop;
  cout << "The population is composed of " << npop << " members" << endl;

  ReadInput >> ncities;
  ReadInput >> generations;
  ReadInput >> SorC;
  ReadInput >> l;

  cout << "There are " << ncities << " cities randomly placed";
  if(SorC==0) cout << " inside a square of side " << l << endl << endl;
  else cout << " on a circumference of radius " << l << endl << endl;


  ReadInput.close();

}

void Initialize()
{

  //creation of the map
  cities = new Map(ncities, SorC, rnd, l);


  //initialization of the probabilities
  p_mutation = new double[mut];
  for (int i=0; i<mut; i++) p_mutation[i]=.1;//rnd->Rannyu(0.,0.1);

  p_crossover=.6;//rnd->Rannyu(0.5,1.);

  //initialization of the starting population
  for (int i=0; i<npop; i++){
    Chromosome x(ncities,mut,p_crossover,p_mutation,rnd);
    population.push_back(x);
    newgen.push_back(x);
    if(x.CheckBonds()!=0) return;
    }

  //print the map
  ofstream WriteMap;
  WriteMap.open("Map.txt");
  for (int i=0; i<ncities; i++){
    WriteMap << i << "\t" << cities->GetX(i) << "\t" << cities->GetY(i) << endl;
  }
  WriteMap.close();

}

int Selection(){

  //occhio: non è necessario l'ordinamento!!

  if(is_sorted(population.begin(),population.end(),order) == false)
    sort(population.begin(),population.end(),order);


  //the probability associated to each chromosome is p ~ 1/(fitness^2)
  vector<double> cdf;   //discrete cumulative density function
  double counter=0;

  for (int i=0; i<npop; i++){
    //double f = exp(-population.at(i).Fitness_L2(*cities));
    double f = 1./pow(population.at(i).Fitness_L2(*cities),2.);
    //p[i]=f;
    counter += f;
    cdf.push_back(counter);
  }

  // for (int i=0; i<npop; i++) p[i] /= counter;

  for (int i=0; i<npop; i++){ cdf.at(i) /= (double)counter;} // cout << i << " " << cdf.at(i) << endl; }

  double r = rnd->Rannyu();
  //cout << r << endl;

  for (int i=0; i<npop; i++)
    if(r<cdf.at(i)) return i;

  return -1;
}

void Crossover(Chromosome father, Chromosome mother, Chromosome* son, Chromosome* daughter ){

  //Achtung!! son and daughter must point Chromosome objects equal to father and mother!!

  double rr = rnd->Rannyu();
  if (rr>p_crossover){
    return;
  }

  int j1=0, j2=0;

  double r = rnd->Rannyu(0,ncities);

  int cut = floor(r);     //decide where to cut the chromosome

  for (int i=cut; i<ncities; i++){    //fill the rest of son with the genetic heritage of the mother
    while (j1 != ncities){
      if( find(son->begin(), son->begin()+cut, mother.GetValue(j1)) == son->begin()+cut ){
        //if j-th city of the mother is not in the sequence of the son
        son->SetValue(i,mother.GetValue(j1));
        //cout << "Eccolo! " << i << " " << j1 << " " << son->GetValue(i) << endl;
        j1++;
        break;
      }
      j1++;
    }
  }

  for (int i=cut; i<ncities; i++){    //fill the rest of daughter with the genetic heritage of the father
    while (j2 != ncities){
      if( find(daughter->begin(), daughter->begin()+cut, father.GetValue(j2)) == daughter->begin()+cut ){
        //if j-th city of the father is not in the sequence of the daughter
        daughter->SetValue(i,father.GetValue(j2));
        j2++;
        break;
      }
      j2++;
    }
  }

}


bool order(Chromosome a, Chromosome b){
  double f1 = a.Fitness_L2(*cities);
  double f2 = b.Fitness_L2(*cities);
  return (f1<f2);
}

void Deallocate(){
  delete cities;
  delete[] p_mutation;
}
