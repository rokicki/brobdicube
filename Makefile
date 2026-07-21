UNAME := $(shell uname)

CXX      := g++
CXXFLAGS := -march=native -Wall -Wsign-compare -g -std=c++20 -O3

ifeq ($(UNAME), Linux)
  LDFLAGS := -lpthread
else
  LDFLAGS :=
endif

.PHONY: all clean

all: brobdicube pareto findpermnet algutil

brobdicube: brobdicube.cpp fpng.cpp
	$(CXX) $(CXXFLAGS) -o $@ brobdicube.cpp -w fpng.cpp $(LDFLAGS)

pareto: pareto.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

findpermnet: findpermnet.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

algutil: algutil.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f brobdicube pareto findpermnet algutil
	rm -rf test
ifeq ($(UNAME), Darwin)
	rm -rf brobdicube.dSYM pareto.dSYM findpermnet.dSYM algutil.dSYM
endif
