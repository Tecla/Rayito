
all: rayito

rayito: main.o
	g++ -o rayito main.o

main.o: main.cpp rayito.h
	g++ -c main.cpp -o main.o -O3 -Wall

clean:
	rm -f main.o rayito out.ppm out.pfm

