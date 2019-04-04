CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: base_sim.exe esercizio4.1.exe

base_sim.exe : MolDyn_NVE.o random.o 
	$(CC) MolDyn_NVE.o random.o -o base_sim.exe
esercizio4.1.exe : esercizio4.1.o random.o 
	$(CC) esercizio4.1.o random.o  -o esercizio4.1.exe

MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
esercizio4.1.o : esercizio4.1.cpp MolDyn_NVE.h
	$(CC) -c esercizio4.1.cpp -o esercizio4.1.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.exe
