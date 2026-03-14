SHELL := /bin/bash

CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp
LDFLAGS = -L/usr/local/lib -lgsl -lgslcblas
INCLUDE = -I../../src

SRC = $(wildcard *.cpp) \
      ../../src/gsl_wrapper/Integration1DimGSL.cpp \
      ../../src/integration_methods/Integration1DimNewtonCotes.cpp

OBJ = $(SRC:.cpp=.o)
DEP = $(OBJ:.o=.d)

TARGET = test.out

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(INCLUDE) $(LDFLAGS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

-include $(DEP)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)
