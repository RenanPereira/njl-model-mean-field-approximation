SHELL := /bin/bash
CXX = g++ -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp

INCLUDE_DIRS = -Isrc

DEPS = Integration1DimNewtonCotes.h

OBJ =  tests.o \
	   Integration1DimNewtonCotes.o


%.o: %.cpp $(DEPS)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

run: $(OBJ) 
	$(CXX) $(OBJ) $(INCLUDE_DIRS)  -o tests.out

clean: 
	rm -f tests.out $(OBJ)
