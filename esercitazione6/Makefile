CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio6.exe : esercizio6.o random.o
	$(CC) random.o esercizio6.o -o esercizio6.exe


esercizio6.o : esercizio6.cpp 
	$(CC) -c esercizio6.cpp -o esercizio6.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)


clean :
	rm *.o esercizio6.exe
