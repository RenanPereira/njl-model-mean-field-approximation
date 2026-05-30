#include <iostream>
#include "TestIniFileParser.h"

int main() 
{
    bool allTestsPassed = true;
    
    // Run each test and collect the result
    allTestsPassed &= testGetValue();
    allTestsPassed &= testGetInt();
    allTestsPassed &= testGetDouble(1E-10);
    allTestsPassed &= testGetSectionData();
    allTestsPassed &= testGetSections();

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
