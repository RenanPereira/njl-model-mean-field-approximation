SHELL := /bin/bash
CXX = g++ -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp

INCLUDE_DIRS = -I. -I../integration_methods

DEPS = Integration1DimGSL.h \
	   ../integration_methods/Integration1DimNewtonCotes.h

OBJ =  tests.o \
	   Integration1DimGSL.o \
	   ../integration_methods/Integration1DimNewtonCotes.o


%.o: %.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

%.o: integration_methods/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

run: $(OBJ) 
	$(CXX) $(OBJ) $(INCLUDE_DIRS)  -L/usr/local/lib -lgsl -lgslcblas -o tests.out

clean: 
	rm -f tests.out $(OBJ)
