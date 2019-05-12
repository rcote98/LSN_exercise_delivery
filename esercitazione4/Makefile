CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio4.exe
	
esercizio4.exe : esercizio4.o random.o 
	$(CC) esercizio4.o random.o  -o esercizio4.exe

esercizio4.o : esercizio4.cpp esercizio4.h
	$(CC) -c esercizio4.cpp -o esercizio4.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.exe
