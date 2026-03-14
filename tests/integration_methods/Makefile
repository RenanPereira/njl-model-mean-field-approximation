SHELL := /bin/bash

CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++17 -fopenmp
LDFLAGS = -fopenmp -L/usr/local/lib
LDLIBS = -lgsl -lgslcblas
INCLUDE = -I../../src

SRC = $(wildcard *.cpp) \
      ../../src/integration_methods/Integration1DimNewtonCotes.cpp

OBJ = $(SRC:.cpp=.o)
DEP = $(OBJ:.o=.d)

TARGET = test.out

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

-include $(DEP)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)
