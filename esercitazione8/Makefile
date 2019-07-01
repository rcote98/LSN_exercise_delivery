CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio8.exe

esercizio8.exe : esercizio8.o random.o metropolis.o sampler.o
	$(CC) esercizio8.o random.o metropolis.o sampler.o -o esercizio8.exe

esercizio8.o : esercizio8.cpp
	$(CC) -c esercizio8.cpp -o esercizio8.o $(CFLAGS)
	
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)

sampler.o : sampler.cpp sampler.h
	$(CC) -c sampler.cpp -o sampler.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.exe
