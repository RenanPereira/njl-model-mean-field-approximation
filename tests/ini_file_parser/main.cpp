#include <iostream>
#include "TestIniFileParser.h"

int main() 
{   
    bool allTestsPassed = TestIniFileParser::runAllTests();

    if (allTestsPassed) 
    {
        std::cout << "\nAll tests passed!\n";
        return 0;
    } 
    else 
    {
        std::cout << "\nSome tests failed!\n";
        return 1;
    }
}
