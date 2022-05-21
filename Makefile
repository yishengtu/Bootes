CC = g++
CFLAGS = -lhdf5 -lhdf5_cpp

compile_main.cpp: src/main.cpp
	$(CC) src/main.cpp $(CFLAGS) -o bin/bootes.out -O3 -fopenmp

