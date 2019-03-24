CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio2.1.exe esercizio2.2.exe

esercizio2.1.exe : esercizio2.1.o random.o
	$(CC) random.o esercizio2.1.o -o esercizio2.1.exe
esercizio2.2.exe : esercizio2.2.o random.o
	$(CC) random.o esercizio2.2.o -o esercizio2.2.exe
esercizio1.1.o : esercizio2.1.cpp
	$(CC) -c esercizio2.1.cpp -o esercizio2.1.o $(CFLAGS)
esercizio2.2.o : esercizio2.2.cpp
	$(CC) -c esercizio2.2.cpp -o esercizio2.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio2.1.exe esercizio2.2.exe seed.out
