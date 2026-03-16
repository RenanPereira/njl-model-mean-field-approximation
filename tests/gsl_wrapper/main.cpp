#include <iostream>
#include "TestIntegration1DimGSL.h"

int main() 
{
    bool allTestsPassed = true;

    // Run each test and collect the result
    allTestsPassed &= hardcodedTestIntegration1DimGSL(1E-8);

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
