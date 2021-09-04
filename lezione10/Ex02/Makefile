CC = g++
CFLAGS = -Wall -O3 --std=c++11



main.exe : main.o random.o classes.o
	mpicxx random.o main.o classes.o -o main.exe 
main.o : main.cxx
	mpicxx -c main.cxx -o main.o $(CFLAGS)
random.o : random.cpp random.h
	mpicxx -c random.cpp -o random.o $(CFLAGS)
classes.o: classes.cxx classes.h random.h
	mpicxx -c -o classes.o classes.cxx $(CFLAGS)


clean :
	rm *.o main.exe seed.out best_half.dat best_path.dat


