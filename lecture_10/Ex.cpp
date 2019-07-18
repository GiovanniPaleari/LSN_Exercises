#include "mpi.h" 
#include <iostream> 
using namespace std;
int main(int argc, char* argv[])
{
MPI_Init(&argc,&argv);
int size, myrank;
MPI_Comm_size(MPI_COMM_WORLD,&size);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
int my_values[3];
for(int i=0;i<3;i++) if(myrank==0) my_values[i]=i+1; else my_values[i]=0;
cout<< "Prima: "<< my_values[0]<< " "<< my_values[1]<< " "<< my_values[2]<< " per il processo "<< myrank<< endl; MPI_Bcast(my_values,3,MPI_INTEGER,0,MPI_COMM_WORLD); cout<< "Dopo: "<< my_values[0]<< " "<< my_values[1]<< " "<< my_values[2]<< " per il processo "<< myrank<< endl;
MPI_Finalize();
return 0; }
