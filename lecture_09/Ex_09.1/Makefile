CC = g++
CFLAGS = -Wall -O3 --std=c++11

GA_tsp.exe : GA_tsp.o random.o Chromosome.o Map.o
	$(CC) random.o Chromosome.o Map.o GA_tsp.o -o GA_tsp.exe
GA_tsp.o : GA_tsp.cpp
	$(CC) -c GA_tsp.cpp -o GA_tsp.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Chromosome.o : Chromosome.cpp Chromosome.h GA_tsp.h
	$(CC) -c Chromosome.cpp -o Chromosome.o $(CFLAGS)
Map.o : Map.cpp Map.h
	$(CC) -c Map.cpp -o Map.o $(CFLAGS)
clean :
	rm *.o GA_tsp.exe
