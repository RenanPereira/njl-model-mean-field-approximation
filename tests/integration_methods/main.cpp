#include <iostream>

#include "TestIntegration1DimNewtonCotes.h"

int main() 
{
    bool allTestsPassed = true;

    // Run each test and collect the result
    allTestsPassed &= testIntegration1DimNewtonCotes(1E-3);

    if (allTestsPassed) 
    {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } 
    else 
    {
        std::cout << "Some tests failed!" << std::endl;
        return 1;
    }
}