CXX = g++
CXXFLAGS = -O2 -std=c++11 -g -Wall

all: align sapling_example binarysearch

align: align.o sapling_api.h ssw_cpp.o ssw.o
	$(CXX) $(CXXFLAGS) -o align align.o ssw_cpp.o ssw.o

clean:
	rm -f align sapling_example binarysearch *.o
