SHELL := /bin/bash

CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++17 -fopenmp
LDFLAGS = -fopenmp -L/usr/local/lib
LDLIBS = -lgsl -lgslcblas
INCLUDE_DIRS = -Isrc

SRC := $(shell find src -name '*.cpp')
OBJ := $(patsubst src/%.cpp,obj/%.o,$(SRC))
DEP := $(OBJ:.o=.d)

TARGET = bin/nambuJonaLasinioModel.out

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -MMD -MP -c $< -o $@

-include $(DEP)

.PHONY: clean

clean:
	rm -f $(TARGET) $(OBJ) $(DEP)
