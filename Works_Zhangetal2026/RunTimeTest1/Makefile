# Makefile for the TLI runtime benchmark on the MPI cluster.

SHELL := /bin/bash

CXX = mpicxx

STD     = -std=c++17
OPT     = -O2 -m64
INC     = -I.

CXXFLAGS = $(OPT) $(STD)  $(INC)
LDFLAGS  = 

SRC  = Test_RuntimeTest.cpp
PROG = Test_RuntimeTest

HEADERS = \
  _integral.h \
  _TLI_integralkernel.h \
  _TLI_estimate.h

.PHONY: all clean cleanObj

all: $(PROG)

$(PROG): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) "$(SRC)" -o $@ $(LDFLAGS)

clean:
	rm -f $(PROG) *.o *.obj *.pdb *.dwf *.exe

cleanObj:
	rm -f *.o *.obj
