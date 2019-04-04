CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio3.1.exe

esercizio3.1.exe : esercizio3.1.o random.o
	$(CC) random.o esercizio3.1.o -o esercizio3.1.exe
esercizio1.1.o : esercizio3.1.cpp
	$(CC) -c esercizio3.1.cpp -o esercizio3.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio3.1.exe seed.out
