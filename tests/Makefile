# Cross-checks src/dmft/* against the original legacy implementations on
# identical inputs (see test_linalg.cpp for why: the original program's RNG
# is seeded from time(NULL), so two runs of the unmodified binary already
# disagree with each other -- a whole-binary before/after diff isn't a
# meaningful check here).
CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -fopenmp
LEGACY_CXXFLAGS = -O2 -fopenmp -std=gnu++14  # see legacy/Makefile

all: test_linalg test_impurity_raw

test_linalg: test_linalg.o legacy_bridge.o ../src/dmft/linalg.o
	$(CXX) -fopenmp -o $@ $^

test_linalg.o: test_linalg.cpp ../src/dmft/linalg.hpp legacy_bridge.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -c $< -o $@

legacy_bridge.o: legacy_bridge.cpp legacy_bridge.hpp ../legacy/headers.new.h
	$(CXX) $(LEGACY_CXXFLAGS) -c $< -o $@

test_impurity_raw: test_impurity_raw.o impurity_raw_bridge.o ../src/dmft/impurity_raw.o ../src/dmft/linalg.o
	$(CXX) -fopenmp -o $@ $^

test_impurity_raw.o: test_impurity_raw.cpp ../src/dmft/impurity_raw.hpp impurity_raw_bridge.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -c $< -o $@

impurity_raw_bridge.o: impurity_raw_bridge.cpp impurity_raw_bridge.hpp ../legacy/headers.new.h ../legacy/parameters.cpp ../legacy/solver.cpp
	$(CXX) $(LEGACY_CXXFLAGS) -c $< -o $@

../src/dmft/linalg.o: ../src/dmft/linalg.cpp ../src/dmft/linalg.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -c $< -o $@

../src/dmft/impurity_raw.o: ../src/dmft/impurity_raw.cpp ../src/dmft/impurity_raw.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -c $< -o $@

check: all
	./test_linalg
	./test_impurity_raw

clean:
	rm -f test_linalg test_impurity_raw *.o ../src/dmft/*.o
