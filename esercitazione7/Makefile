CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio7.exe : esercizio7.o random.o
	$(CC) random.o esercizio7.o -o esercizio7.exe
esercizio7.o : esercizio7.cpp
	$(CC) -c esercizio7.cpp -o esercizio7.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio7.exe
