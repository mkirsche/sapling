CXX = g++
CXXFLAGS = -O2 -std=c++11 -g -Wall

sapling: align.o sapling_api.h ssw_cpp.o ssw.o
	$(CXX) $(CXXFLAGS) -o align align.o sapling_api.h ssw_cpp.o ssw.o

clean:
	rm align *.o
