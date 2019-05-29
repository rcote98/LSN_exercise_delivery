CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio8.exe

esercizio8.exe : esercizio8.o random.o
	$(CC) esercizio8.o random.o -o esercizio8.exe

esercizio8.o : esercizio8.cpp esercizio8.h
	$(CC) -c esercizio8.cpp -o esercizio8.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.exe *.csv
