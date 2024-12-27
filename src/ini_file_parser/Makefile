SHELL := /bin/bash
CXX = g++ -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp

INCLUDE_DIRS = -Isrc

DEPS = IniFileParser.h

OBJ =  main_tests.o \
	   IniFileParser.o


%.o: %.cpp $(DEPS)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

run: $(OBJ) 
	$(CXX) $(OBJ) $(INCLUDE_DIRS)  -o run_tests.out

clean: 
	rm -f run_tests.out $(OBJ)
