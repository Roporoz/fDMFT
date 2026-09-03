# Reference build for the original (unmodified) solver.
# Requires -std=gnu++14: with C++17+, the global `double beta` in
# parameters.cpp collides with std::beta() pulled in by <cmath>.
CXX = g++
CXXFLAGS = -O2 -std=gnu++14 -fopenmp

dmft_legacy: main.cpp headers.new.h parameters.cpp solver.cpp
	$(CXX) $(CXXFLAGS) -o dmft_legacy main.cpp

clean:
	rm -f dmft_legacy *.dat *.log
