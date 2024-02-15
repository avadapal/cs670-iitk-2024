CXX=g++
CXXFLAGS=--std=c++17 -march=native -O3 -g -Wno-ignored-attributes -fopenmp -DBOOST_ERROR_CODE_HEADER_ONLY
LIBS=-lbsd -lboost_system -lboost_thread -pthread

all: dpf

dpf: dpf.cpp dpf.h prg.h block.h
	g++ -ggdb -O3 -march=native -std=c++17 -Wall -pedantic -o dpf dpf.cpp -lbsd -Wno-ignored-attributes

clean:
	rm -f  dpf
