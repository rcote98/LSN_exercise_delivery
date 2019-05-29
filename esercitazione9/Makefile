CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio9.exe

esercizio9.exe : esercizio9.o random.o traveling.o
	$(CC) esercizio9.o random.o traveling.o -o esercizio9.exe

esercizio9.o : esercizio9.cpp
	$(CC) -c esercizio9.cpp -o esercizio9.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

traveling.o : traveling.cpp traveling.h
	$(CC) -c traveling.cpp -o traveling.o $(CFLAGS)

clean:
	rm *.o *.exe output.*.dat
