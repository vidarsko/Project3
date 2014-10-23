# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -larmadillo -O1
file = test_trialwavefunction

# ****************************************************
# Targets needed to bring the executable up to date

$(file).x: $(file).cpp project3lib.cpp
	$(CXX)  -o $(file).x  project3lib.cpp $(file).cpp $(CXXFLAGS)