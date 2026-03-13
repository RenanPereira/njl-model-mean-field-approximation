CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp
INCLUDE = -I../../src

SRC = main.cpp \
	  TestIntegration1DimNewtonCotes.cpp \
      ../../src/integration_methods/Integration1DimNewtonCotes.cpp

TARGET = test.out

$(TARGET):
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
