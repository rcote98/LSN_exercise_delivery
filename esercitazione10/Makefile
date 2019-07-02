CC = g++
PC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

all: esercizio10.1.exe esercizio10.2.exe

esercizio10.1.exe : esercizio10.1.o random.o annealing.o
	$(CC) esercizio10.1.o random.o annealing.o -o esercizio10.1.exe

esercizio10.2.exe : esercizio10.2.o random.o parallel.o
	$(PC) esercizio10.2.o random.o parallel.o -o esercizio10.2.exe

esercizio10.1.o : esercizio10.1.cpp
	$(CC) -c esercizio10.1.cpp -o esercizio10.1.o $(CFLAGS)

esercizio10.2.o : esercizio10.2.cpp
	$(PC) -c esercizio10.2.cpp -o esercizio10.2.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

annealing.o : annealing.cpp annealing.h
	$(PC) -c annealing.cpp -o annealing.o $(CFLAGS)

parallel.o : parallel.cpp parallel.h
	$(PC) -c parallel.cpp -o parallel.o $(CFLAGS)

clean:
	rm *.o *.exe
