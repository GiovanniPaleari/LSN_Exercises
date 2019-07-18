#include "SA_tsp.h"

using namespace std;

int main()
{
  Input();
  Initialize();
  best.open("Path.txt");

  cout << "Length of the initial path = " << x->Fitness_L1(*cities) << endl;

  x->PrintAdjList("InitialPath.txt");

  int cc=0;

  for(int i=0; i<nstep.size(); i++){
    for(int j=0; j<nstep.at(i); j++){
      Move(i);
    }
    cout << "------------------------------" << endl;
    cout << "Beta_" << i << " = " << beta.at(i) << endl;
    cout << "Nstep_" << i << " = " << nstep.at(i) << endl;
    cout << "Acceptance rate = " << accepted/attempted << endl;
    cout << "Length of the actual path = " << x->Fitness_L1(*cities) << endl;
    cout << i << "\t" << x->Fitness_L1(*cities) << endl;
    best << i << "\t" << x->Fitness_L1(*cities) << endl;
    accepted=0;
    attempted=0;


  }


  ConfFinal();

  //preventing memory leak
  delete rnd;
  delete cities;
  delete[] p_mutation;
  delete x;

  best.close();
  return 0;
}

void Input(void)
{

  cout << "Simulated Annealing Algorithm to optimize the Traveling Salesman Problem" << endl << endl;
  cout << "A Chromosome is a possible closed path which visit each city only once" << endl << endl;

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

  ifstream ReadInput, ReadAnnSched;
  ofstream WriteAnnSched;

  ReadInput.open("input.dat");

  ReadInput >> ncities;
  ReadInput >> SorC;
  ReadInput >> l;

  cout << "There are " << ncities << " cities randomly placed";
  if(SorC==0) cout << " inside a square of side " << l << endl << endl;
  else cout << " on a circumference of radius " << l << endl << endl;


  ReadInput.close();

  WriteAnnSched.open("AnnealingSchedule.txt");
  double t=1;
  int s=1000;
  double alpha = 1.1;
  for (int i=0; i<100; i++){
    beta.push_back(t*alpha);
    nstep.push_back(s);
    WriteAnnSched << t*alpha << "\t" << s << endl;
    t=t*alpha;
  }


}

void Initialize()
{

  accepted=0;
  attempted=0;

  //creation of the map
  cities = new Map(ncities, SorC, rnd, l);


  //initialization of the probabilities (even if they are useless here)
  p_mutation = new double[mut];
  for (int i=0; i<mut; i++) p_mutation[i]=.1;//rnd->Rannyu(0.,0.1);

  p_crossover=.6;//rnd->Rannyu(0.5,1.);

  //initialization of the starting chromosome
  x = new Chromosome(ncities,mut,p_crossover,p_mutation,rnd);
  if(x->CheckBonds()!=0) return;


  //print the map
  ofstream WriteMap;
  WriteMap.open("Map.txt");
  for (int i=0; i<ncities; i++){
    WriteMap << i << "\t" << cities->GetX(i) << "\t" << cities->GetY(i) << endl;
  }
  WriteMap.close();

//useless for a Simulated Annealing simulation
  for (int i=0; i<npop; i++){
    Chromosome y(ncities,mut,p_crossover,p_mutation,rnd);
    population.push_back(y);
    newgen.push_back(y);
    if(y.CheckBonds()!=0) return;
    }

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

void Move(int sched)
{
  x->Shift(); //in order to make the algorithm ergodic

  int n = floor(rnd->Rannyu(0.,ncities));
  int m=n;
  while (m==n) m = floor(rnd->Rannyu(0.,ncities));
  //int m = delta;

  Chromosome tmp(x);

  reverse(tmp.begin()+min(n,m), tmp.begin()+max(n,m));    //mossa splendida, capire meglio perchè


  double p = exp(-beta.at(sched)*(tmp.Fitness_L1(*cities) - x->Fitness_L1(*cities)));
  double r = rnd->Rannyu();

  if (r<p){
    *x=tmp;
    accepted = accepted + 1.0;
    //cout << "Accepted: swap between " << n << " and " << m << endl;
  }
  attempted = attempted + 1.0;

}

void ConfFinal(void)
{
  cout << "Print adjacency list to file BestPath.txt " << endl << endl;
  x->PrintAdjList("BestPath.txt");

  rnd->SaveSeed();
}
