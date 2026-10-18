SHELL := /bin/bash

CXX = mpicxx
STD = -std=c++17
OPT = -O2 -m64
CPPFLAGS = -I.
CXXFLAGS = $(OPT) $(STD)  $(CPPFLAGS)
LDFLAGS = 

PROGS = \
  Work1_GravityEstimateMPI \
  Work1_MagneticEstimateMPI

HEADERS = \
  _TLI_estimate.h \
  _TLI_integralkernel.h \
  _integral.h

.PHONY: all clean

all: $(PROGS)

Work1_GravityEstimateMPI: Work1_GravityEstimateMPI.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS)

Work1_MagneticEstimateMPI: Work1_MagneticEstimateMPI.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS)

clean:
	rm -f $(PROGS)
	rm -f *.o *.obj *.pdb *.dwf *.exe
