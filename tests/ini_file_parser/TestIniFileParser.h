#ifndef TESTINIFILEPARSER_H
#define TESTINIFILEPARSER_H

class TestIniFileParser
{
private:
    static bool check(bool , const std::string& );

    static void printTestResult(bool , const std::string& );

    static bool testGetValue();
    static bool testGetInt();
    static bool testGetDouble(double );
    static bool testGetSectionData();
    static bool testGetSections();
    static bool testValidatePositiveDouble();
    static bool testIsKeyPresent();
    static bool testValidatePositiveInteger();
    static bool testValidateNonNegativeDouble();
    static bool testValidateRequiredSections();
    static bool testValidateRequiredKeys();
    static bool testGetBool();
    static bool testValidateBool();

public:
    static bool runAllTests();
};


#endif