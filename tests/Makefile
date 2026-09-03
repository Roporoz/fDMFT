# Cross-checks src/dmft/linalg.* against the original legacy/headers.new.h
# implementations on identical inputs (see test_linalg.cpp for why: the
# original program's RNG is seeded from time(NULL), so two runs of the
# unmodified binary already disagree with each other -- a whole-binary
# before/after diff isn't a meaningful check here).
CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -fopenmp
LEGACY_CXXFLAGS = -O2 -fopenmp -std=gnu++14  # see legacy/Makefile

test_linalg: test_linalg.o legacy_bridge.o ../src/dmft/linalg.o
	$(CXX) -fopenmp -o $@ $^

test_linalg.o: test_linalg.cpp ../src/dmft/linalg.hpp legacy_bridge.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -c $< -o $@

legacy_bridge.o: legacy_bridge.cpp legacy_bridge.hpp ../legacy/headers.new.h
	$(CXX) $(LEGACY_CXXFLAGS) -c $< -o $@

../src/dmft/linalg.o: ../src/dmft/linalg.cpp ../src/dmft/linalg.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -c $< -o $@

check: test_linalg
	./test_linalg

clean:
	rm -f test_linalg *.o ../src/dmft/*.o
