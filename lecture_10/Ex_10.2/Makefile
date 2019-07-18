CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

SA_tsp.exe : SA_tsp.o random.o Chromosome.o Map.o
	$(CC) random.o Chromosome.o Map.o SA_tsp.o -o SA_tsp.exe
SA_tsp.o : SA_tsp.cpp
	$(CC) -c SA_tsp.cpp -o SA_tsp.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Chromosome.o : Chromosome.cpp Chromosome.h SA_tsp.h
	$(CC) -c Chromosome.cpp -o Chromosome.o $(CFLAGS)
Map.o : Map.cpp Map.h
	$(CC) -c Map.cpp -o Map.o $(CFLAGS)
clean :
	rm *.o SA_tsp.exe

#exe: mpiexec -np 4 --oversubscribe SA_tsp.exe 
