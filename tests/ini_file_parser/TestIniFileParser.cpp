#include <iostream>
#include <cmath>
#include "ini_file_parser/IniFileParser.h"
#include "TestIniFileParser.h"

// Check function that simply prints an error message and returns false on failure
bool TestIniFileParser::check(bool condition, const std::string& message) 
{
    if (!condition) 
    {
        std::cerr << "Test failed: " << message << std::endl;
        return false;
    }
    else
    {
        return true;
    }
}

void TestIniFileParser::printTestResult(bool testPassed, const std::string& testFunctionName)
{
    if (testPassed)
    {
        std::cout << testFunctionName  << " passed.\n";
    }
    else
    {
        std::cout << testFunctionName  << " failed.\n";
    }
}

bool TestIniFileParser::testGetValue() 
{
    std::cout << "\nRunning " << __func__ << "...\n";
    
    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[General]\n";
    out << "version = 1.2.3\n";
    out << "name = TestApp\n";
    out << "[Settings]\n";
    out << "fullscreen = true\n";
    out << "width = 1920\n";
    out << "height = 1080\n";
    out.close();

    IniFileParser parser(testIniFile);

    // Check basic value retrieval
    bool passed = true;
    passed &= check(parser.getValue("General", "version") == "1.2.3", "version value mismatch");
    passed &= check(parser.getValue("General", "name") == "TestApp", "name value mismatch");
    passed &= check(parser.getValue("Settings", "fullscreen") == "true", "fullscreen value mismatch");

    // Check non-existing value
    passed &= check(parser.getValue("General", "nonexistent") == "", "expected empty string for nonexistent key");

    // Clean up test file
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);
    
    return passed;
}

bool TestIniFileParser::testGetInt() 
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Display]\n";
    out << "width = 1920\n";
    out << "height = 1080\n";
    out.close();

    IniFileParser parser(testIniFile);

    // Check integer values
    bool passed = true;
    passed &= check(parser.getInt("Display", "width") == 1920, "width value mismatch given string section");
    passed &= check(parser.getInt("Display", "height") == 1080, "height value mismatch given string section");

    // Check non-existing integer (should return 0)
    passed &= check(parser.getInt("Display", "nonexistent") == 0, "expected 0 for nonexistent integer key given string section");
    passed &= check(parser.getInt("Display", "nonexistent", 100) == 100, "expected 100 (default value) for nonexistent integer key given string section");

    // Create a map to simulate a section in the INI file
    std::map<std::string, std::string> testSection = {
        {"width", "1920"},
        {"height", "1080"}
    };

    // Test existing keys
    passed &= check(parser.getInt(testSection, "width") == 1920, "width value mismatch given map");
    passed &= check(parser.getInt(testSection, "height") == 1080, "height value mismatch given map");

    // Test non-existing key
    passed &= check(parser.getInt(testSection, "nonexistent") == 0, "expected 0 for nonexistent integer key given map");
    passed &= check(parser.getInt(testSection, "nonexistent", 10) == 10, "expected 10 (default value) for nonexistent integer key given map");

    // Clean up test file
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testGetDouble(double absError) 
{   
    std::cout << "\nRunning " << __func__ << "...\n";
    
    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Settings]\n";
    out << "volume = 0.75\n";
    out << "brightness = 0.9\n";
    out.close();

    IniFileParser parser(testIniFile);

    // Test existing keys
    bool passed = true;
    passed &= check(fabs(parser.getDouble("Settings", "volume") - 0.75) < absError, "volume value mismatch given string section");
    passed &= check(fabs(parser.getDouble("Settings", "brightness") - 0.9 ) < absError, "brightness value mismatch given string section");

    // Test non-existing keys
    passed &= check(parser.getDouble("Settings", "nonexistent") < absError, "expected 0.0 for nonexistent double key given string section");
    passed &= check(fabs(parser.getDouble("Settings", "nonexistent", 3.14) - 3.14) < absError, "expected 3.14 (default value) for nonexistent double key given string section");

    // Create a map to simulate a section in the INI file
    std::map<std::string, std::string> testSection = {
        {"volume", "0.75"},
        {"brightness", "0.9"}
    };

    // Test existing keys
    passed &= check(fabs(parser.getDouble(testSection, "volume") - 0.75 ) < absError, "volume mismatch given map");
    passed &= check(fabs(parser.getDouble(testSection, "brightness") - 0.9 ) < absError, "brightness mismatch given map");

    // Test non-existing key
    passed &= check(fabs(parser.getDouble(testSection, "nonexistent") - 0.0 ) < absError, "expected 0.0 for nonexistent key given map");
    passed &= check(fabs(parser.getDouble(testSection, "nonexistent", 3.14) - 3.14 ) < absError, "expected 3.14 (default value) for nonexistent key given map");

    // Clean up test file
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testGetSectionData() 
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[User1]\n";
    out << "name = Alice\n";
    out << "age = 30\n";
    out << "[User2]\n";
    out << "name = Bob\n";
    out << "age = 25\n";
    out << "[User3]\n";
    out << "name = Carol\n";
    out << "age = 27\n";
    out.close();

    IniFileParser parser(testIniFile);

    // Test getSectionsData() with a common prefix "User"
    auto sections = parser.getSectionsData("User");

    bool passed = true;
    passed &= check(sections.size() == 3, "Expected 3 sections with prefix 'User'");

    // Check values in each section
    passed &= check(sections[0].at("age") == "30", "User1 age mismatch");
    passed &= check(sections[1].at("name") == "Bob", "User2 name mismatch");
    passed &= check(sections[2].at("age") == "27", "User3 age mismatch");

    // Clean up test file
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testGetSections()
{
    std::cout << "\nRunning " << __func__ << "...\n";
    
    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[User1]\n";
    out << "name = Alice\n";
    out << "age = 30\n";
    out << "[User2]\n";
    out << "name = Bob\n";
    out << "age = 25\n";
    out << "[User3]\n";
    out << "name = Carol\n";
    out << "age = 27\n";
    out.close();

    IniFileParser parser(testIniFile);

    // Test getSections() with a common prefix "User"
    auto sections = parser.getSections("User");

    bool passed = true;

    // Check if 3 sections are returned
    passed &= check(sections.size() == 3, "Expected 3 sections with prefix 'User'");

    // Check section names and corresponding values
    passed &= check(sections[0].first == "User1", "User1 section name mismatch");
    passed &= check(sections[0].second.at("name") == "Alice", "User1 name mismatch");
    passed &= check(sections[0].second.at("age") == "30", "User1 age mismatch");

    passed &= check(sections[1].first == "User2", "User2 section name mismatch");
    passed &= check(sections[1].second.at("name") == "Bob", "User2 name mismatch");
    passed &= check(sections[1].second.at("age") == "25", "User2 age mismatch");

    passed &= check(sections[2].first == "User3", "User3 section name mismatch");
    passed &= check(sections[2].second.at("name") == "Carol", "User3 name mismatch");
    passed &= check(sections[2].second.at("age") == "27", "User3 age mismatch");

    // Clean up the test file
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testValidatePositiveDouble()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Settings]\n";
    out << "positiveValue = 0.75\n";
    out << "zeroValue = 0.0\n";
    out << "negativeValue = -5.0\n";
    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    // key exists and value > 0
    passed &= check(
        parser.validatePositiveDouble(
            "Settings",
            "positiveValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveDouble to return true for positive double value"
    );

    // key exists and value == 0
    passed &= check(
        !parser.validatePositiveDouble(
            "Settings",
            "zeroValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveDouble to return false for zero value"
    );

    // key exists and value < 0
    passed &= check(
        !parser.validatePositiveDouble(
            "Settings",
            "negativeValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveDouble to return false for negative value"
    );

    // key does not exist
    passed &= check(
        !parser.validatePositiveDouble(
            "Settings",
            "missingValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveDouble to return false for missing key"
    );
    
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testIsKeyPresent()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Section]\n";
    out << "key = value\n";
    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    // key exist
    passed &= check(
        parser.isKeyPresent(
            "Section",
            "key"),
        "Expected isKeyPresent to return true for key existing in file"
    );

    // key does not exist
    passed &= check(
        !parser.isKeyPresent(
            "Section",
            "missingKey"),
        "Expected isKeyPresent to return false for key not existing in file"
    );

    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testValidatePositiveInteger()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Settings]\n";
    out << "positiveValue = 2\n";
    out << "zeroValue = 0\n";
    out << "negativeValue = -5\n";
    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    // key exists and value > 0
    passed &= check(
        parser.validatePositiveInteger(
            "Settings",
            "positiveValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveInteger to return true for positive integer value"
    );

    // key exists and value == 0
    passed &= check(
        !parser.validatePositiveInteger(
            "Settings",
            "zeroValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveInteger to return false for zero value"
    );

    // key exists and value < 0
    passed &= check(
        !parser.validatePositiveInteger(
            "Settings",
            "negativeValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveInteger to return false for negative integer value"
    );

    // key does not exist
    passed &= check(
        !parser.validatePositiveInteger(
            "Settings",
            "missingValue",
            "Invalid file",
            "Value must be positive"),
        "Expected validatePositiveInteger to return false for missing key"
    );
    
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testValidateNonNegativeDouble()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Settings]\n";
    out << "positiveValue = 0.75\n";
    out << "zeroValue = 0.0\n";
    out << "negativeValue = -5.0\n";
    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    // key exists and value >= 0
    passed &= check(
        parser.validateNonNegativeDouble(
            "Settings",
            "positiveValue",
            "Invalid file",
            "Value must be non negative"),
        "Expected validateNonNegativeDouble to return true for positive double value"
    );

    // key exists and value == 0
    passed &= check(
        parser.validateNonNegativeDouble(
            "Settings",
            "zeroValue",
            "Invalid file",
            "Value must be non negative"),
        "Expected validateNonNegativeDouble to return true for zero value"
    );

    // key exists and value < 0
    passed &= check(
        !parser.validateNonNegativeDouble(
            "Settings",
            "negativeValue",
            "Invalid file",
            "Value must be non negative"),
        "Expected validateNonNegativeDouble to return false for negative value"
    );

    // key does not exist
    passed &= check(
        !parser.validateNonNegativeDouble(
            "Settings",
            "missingValue",
            "Invalid file",
            "Value must be non negative"
        ),
        "Expected validateNonNegativeDouble to return false for missing key"
    );
    
    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testValidateRequiredSections()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Section A]\n";
    out << "keyA1 = valueA1\n";
    out << "[Section B]\n";
    out << "keyB1 = valueB1\n";
    out << "[Section C]\n";
    out << "keyC1 = valueC1\n";
    out << "[Section F]\n";
    out << "[Section G]\n";
    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    // all sections present
    passed &= check(
        parser.validateRequiredSections(
            {"Section A", "Section B", "Section C"}
        ),
        "Expected validateRequiredSections to return true when all required sections are present"
    );

    passed &= check(
        !parser.validateRequiredSections(
            {"Section A", "Section B", "Section C", "Section D", "Section E"}
        ),
        "Expected validateRequiredSections to return false when some required sections are missing"
    );

    passed &= check(
        parser.validateRequiredSections(
            {"Section F", "Section G"}
        ),
        "Expected validateRequiredSections to return true when required sections are present but keys and values are missing"
    );

    passed &= check(
        parser.validateRequiredSections(
            {}
        ),
        "Expected validateRequiredSections to return true when no required sections are present"
    );

    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testValidateRequiredKeys()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Section A]\n";
    out << "keyA1 = valueA1\n";
    out << "keyA2 = valueA2\n";
    out << "[Section B]\n";
    out << "keyB1 = valueB1\n";
    out << "[Section C]\n";
    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    // all required keys present
    passed &= check(
        parser.validateRequiredKeys(
            "Section A",
            {"keyA1", "keyA2"}
        ),
        "Expected validateRequiredKeys to return true when all required keys are present"
    );

    // multiple required keys missing
    passed &= check(
        !parser.validateRequiredKeys(
            "Section B",
            {"keyB2", "keyB3"}
        ),
        "Expected validateRequiredKeys to return false when multiple required keys are missing"
    );

    // empty required key list
    passed &= check(
        parser.validateRequiredKeys(
            "Section C",
            {}
        ),
        "Expected validateRequiredKeys to return true when no required keys are specified"
    );

    // empty required key list
    passed &= check(
        parser.validateRequiredKeys(
            "Section A",
            {}
        ),
        "Expected validateRequiredKeys to return true when no required keys are specified"
    );

    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testGetBool()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Section]\n";

    out << "trueValue = true\n";
    out << "oneValue = 1\n";
    out << "yesValue = yes\n";
    out << "onValue = on\n";

    out << "falseValue = false\n";
    out << "zeroValue = 0\n";
    out << "noValue = no\n";
    out << "offValue = off\n";

    out << "trimmedTrue =   true   \n";
    out << "invalidValue = maybe\n";
    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    // true values
    passed &= check(
        parser.getBool("Section", "trueValue", false) == true,
        "Expected getBool to return true for value 'true'"
    );

    passed &= check(
        parser.getBool("Section", "oneValue", false) == true,
        "Expected getBool to return true for value '1'"
    );

    passed &= check(
        parser.getBool("Section", "yesValue", false) == true,
        "Expected getBool to return true for value 'yes'"
    );

    passed &= check(
        parser.getBool("Section", "onValue", false) == true,
        "Expected getBool to return true for value 'on'"
    );

    // false values
    passed &= check(
        parser.getBool("Section", "falseValue", true) == false,
        "Expected getBool to return false for value 'false'"
    );

    passed &= check(
        parser.getBool("Section", "zeroValue", true) == false,
        "Expected getBool to return false for value '0'"
    );

    passed &= check(
        parser.getBool("Section", "noValue", true) == false,
        "Expected getBool to return false for value 'no'"
    );

    passed &= check(
        parser.getBool("Section", "offValue", true) == false,
        "Expected getBool to return false for value 'off'"
    );

    // trimming
    passed &= check(
        parser.getBool("Section", "trimmedTrue", false) == true,
        "Expected getBool to trim whitespace before parsing"
    );

    // invalid value uses default
    passed &= check(
        parser.getBool("Section", "invalidValue", true) == true,
        "Expected getBool to return default value for invalid boolean string"
    );

    passed &= check(
        parser.getBool("Section", "invalidValue", false) == false,
        "Expected getBool to return supplied default value for invalid boolean string"
    );

    // missing key uses default
    passed &= check(
        parser.getBool("Section", "missingKey", true) == true,
        "Expected getBool to return default value when key is missing"
    );

    passed &= check(
        parser.getBool("Section", "missingKey", false) == false,
        "Expected getBool to return supplied default value when key is missing"
    );

    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::testValidateBool()
{
    std::cout << "\nRunning " << __func__ << "...\n";

    const std::string testIniFile = std::string(__func__) + ".ini";

    std::ofstream out(testIniFile);
    out << "[Section]\n";

    out << "trueValue = true\n";
    out << "oneValue = 1\n";
    out << "yesValue = yes\n";
    out << "onValue = on\n";

    out << "falseValue = false\n";
    out << "zeroValue = 0\n";
    out << "noValue = no\n";
    out << "offValue = off\n";

    out << "invalidValue = maybe\n";
    out << "emptyValue = \n";

    out << "trueWithExtraSpaces =      true     \n";

    out.close();

    IniFileParser parser(testIniFile);

    bool passed = true;

    std::vector<std::pair<std::string, std::string> > validCases =
    {
        {"trueValue",  "true"},
        {"oneValue",   "1"},
        {"yesValue",   "yes"},
        {"onValue",    "on"},
        {"falseValue", "false"},
        {"zeroValue",  "0"},
        {"noValue",    "no"},
        {"offValue",   "off"}
    };

    const std::string conditionMessage = "Value must be 'true', '1', 'yes', 'on', 'false', '0', 'no' or 'off'";

    for (int i = 0; i < int(validCases.size()); ++i)
    {
        passed &= check(
            parser.validateBool(
                "Section",
                validCases[i].first,
                "Invalid file",
                conditionMessage
            ),
            "Expected validateBool to return true for '" + validCases[i].second + "'"
        );
    }
    
    // invalid values
    passed &= check(
        !parser.validateBool(
            "Section", 
            "invalidValue", 
            "Invalid file", 
            conditionMessage
        ),
        "Expected validateBool to return false for invalid boolean string"
    );

    // empty values
    passed &= check(
        !parser.validateBool(
            "Section", 
            "emptyValue", 
            "Invalid file", 
            conditionMessage
        ),
        "Expected validateBool to return false for empty string"
    );

    // missing key
    passed &= check(
        !parser.validateBool(
            "Section", 
            "missingKey", 
            "Invalid file", 
            conditionMessage
        ),
        "Expected validateBool to return false for missing key"
    );

    // trimmed value
    passed &= check(
        parser.validateBool(
            "Section",
            "trueWithExtraSpaces",
            "Invalid file",
            conditionMessage
        ),
        "Expected validateBool to ignore surrounding whitespace"
    );

    std::remove(testIniFile.c_str());

    printTestResult(passed, __func__);

    return passed;
}

bool TestIniFileParser::runAllTests()
{
    bool passed = true;

    passed &= testGetValue();
    passed &= testGetInt();
    passed &= testGetDouble(1e-9);
    passed &= testGetSectionData();
    passed &= testGetSections();
    passed &= testValidatePositiveDouble();
    passed &= testValidatePositiveInteger();
    passed &= testValidateNonNegativeDouble();
    passed &= testValidateRequiredSections();
    passed &= testValidateRequiredKeys();
    passed &= testGetBool();
    passed &= testValidateBool();

    return passed;
}
