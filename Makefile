UNAME := $(shell uname)

CXX      := g++
CXXFLAGS := -march=native -Wall -Wsign-compare -g -std=c++20 -O3 -Isrc

ifeq ($(UNAME), Linux)
  LDFLAGS := -lpthread
else
  LDFLAGS :=
endif

S := src

.PHONY: all clean

all: brobdicube pareto findpermnet algutil test_planmask

brobdicube: $S/brobdicube.cpp $S/xcube.cpp $S/cornersolver.cpp $S/wingsolver.cpp $S/netgen.cpp $S/util.cpp $S/fpng.cpp
	$(CXX) $(CXXFLAGS) -o $@ $S/brobdicube.cpp $S/xcube.cpp $S/cornersolver.cpp $S/wingsolver.cpp $S/netgen.cpp $S/util.cpp -w $S/fpng.cpp $(LDFLAGS)

pareto: $S/pareto.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

findpermnet: $S/findpermnet.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

algutil: $S/algutil.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

test_planmask: $S/test_planmask.cpp $S/planmask.h
	$(CXX) $(CXXFLAGS) -o $@ $S/test_planmask.cpp $(LDFLAGS)

clean:
	rm -f brobdicube pareto findpermnet algutil test_planmask
	rm -rf test
ifeq ($(UNAME), Darwin)
	rm -rf brobdicube.dSYM pareto.dSYM findpermnet.dSYM algutil.dSYM
endif
