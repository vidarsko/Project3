# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -larmadillo -O1
file = test

# ****************************************************
# Targets needed to bring the executable up to date

$(file).x: $(file).cpp project3lib.cpp
	$(CXX) $(CXXFLAGS) -o $(file).x $(file).cpp project3lib.cpp