SHELL := /bin/bash

CXX = mpicxx
STD = -std=c++17
OPT = -O2 -m64
CPPFLAGS = -DTLI_ENABLE_MPI -I.
CXXFLAGS = $(OPT) $(STD)  $(CPPFLAGS)
LDFLAGS = 

PROGS = \
  Work1_MassShellIntegrationMPI \
  Work1_MagnetizationShellIntegrationMPI \
  Work2_MassShellIntegrationVzMPI \
  Work2_MagnetizationShellIntegrationVectorMPI

HEADERS = \
  _tesseroid_MSI.h \
  _TLI_estimate.h \
  _TLI_integralkernel.h \
  _integral.h \
  _integral2.h \
  _sphinterp.h \
  _tesseroid_surfacekernel.h

.PHONY: all clean

all: $(PROGS)

Work1_MassShellIntegrationMPI: Work1_MassShellIntegrationMPI.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS)

Work1_MagnetizationShellIntegrationMPI: Work1_MagnetizationShellIntegrationMPI.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS)

Work2_MassShellIntegrationVzMPI: Work2_MassShellIntegrationVzMPI.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS)

Work2_MagnetizationShellIntegrationVectorMPI: Work2_MagnetizationShellIntegrationVectorMPI.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS)

clean:
	rm -f $(PROGS)
	rm -f *.o *.obj *.pdb *.dwf *.exe
