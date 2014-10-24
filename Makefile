# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -larmadillo -O1
CXXFLAGSCOMP = -larmadillo -O1 -c
file = test_multiparticles

# ****************************************************
# Targets needed to bring the executable up to date

$(file).x: $(file).o project3lib.o
	$(CXX)  -o $(file).x  project3lib.o $(file).o $(CXXFLAGS)

$(file).o: $(file).cpp 
	$(CXX) $(file).cpp $(CXXFLAGSCOMP)

project3lib.o: project3lib.cpp 
	$(CXX) project3lib.cpp $(CXXFLAGSCOMP) 

