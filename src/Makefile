# Builds the rewritten DMFT solver (replaces legacy/Makefile's dmft_legacy).
# Unlike the original, this needs no -std=gnu++14 workaround: there's no
# global named `beta` colliding with std::beta from <cmath>.
CXX = g++
CXXFLAGS = -O2 -std=c++17 -fopenmp

OBJS = dmft/linalg.o dmft/impurity_raw.o dmft/impurity_fine.o dmft/dmft_loop.o main.o

dmft_solver: $(OBJS)
	$(CXX) -fopenmp -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f dmft_solver $(OBJS)
