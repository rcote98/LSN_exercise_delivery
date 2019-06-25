CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio5.exe

esercizio5.exe : esercizio5.o random.o metropolis.o
	$(CC) esercizio5.o random.o metropolis.o -o esercizio5.exe

esercizio5.o : esercizio5.cpp
	$(CC) -c esercizio5.cpp -o esercizio5.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)

clean:
	rm *.o *.exe
