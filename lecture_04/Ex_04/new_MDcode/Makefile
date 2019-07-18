CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : MolDyn_NVE.o random.o
	$(CC)   MolDyn_NVE.o random.o -o main.exe
random.o : random.cpp random.h
		$(CC) -c random.cpp -o random.o $(CFLAGS)
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
