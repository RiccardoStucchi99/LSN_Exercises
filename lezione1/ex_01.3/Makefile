LIBS := `root-config --libs`
INCS := `root-config --cflags`



main.exe : main.o random.o
	g++ random.o main.o -o main.exe ${INCS} ${LIBS}
main.o : main.cxx
	g++ -c main.cxx -o main.o $(CFLAGS) ${INCS} ${LIBS}
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out 
