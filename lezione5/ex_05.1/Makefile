CC = g++
CFLAGS = -Wall -O3 --std=c++11



main.exe : main.o random.o classes.o
	$(CC) random.o main.o classes.o -o main.exe 
main.o : main.cxx
	$(CC) -c main.cxx -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
classes.o: classes.cxx classes.h random.h
	g++ -c -o classes.o classes.cxx $(CFLAGS)


clean :
	rm *.o main.exe seed.out


