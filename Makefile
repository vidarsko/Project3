# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -larmadillo -O2   -fopenmp 
CXXFLAGSCOMP = -larmadillo -O2 -c -fopenmp 
file = test_multiparticles
#file = test_investigate

# ****************************************************
# Targets needed to bring the executable up to date

$(file).x: $(file).o project3lib.o
	$(CXX)  -o $(file).x  project3lib.o $(file).o $(CXXFLAGS)

$(file).o: $(file).cpp project3lib.o
	$(CXX) $(file).cpp project3lib.o $(CXXFLAGSCOMP)

project3lib.o: project3lib.cpp 
	$(CXX) project3lib.cpp $(CXXFLAGSCOMP) 

