#include <iostream>
#include <cmath>
#include "IniFileParser.h"


// Test functions return true if the test passes, false if it fails
bool check(bool , const std::string& );
bool testGetValue();
bool testGetInt();
bool testGetDouble(double );
bool testGetSectionData();
bool testGetSections();


// Check function that simply prints an error message and returns false on failure
bool check(bool condition, const std::string& message) 
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

bool testGetValue() 
{
    const std::string testIniFile = "test.ini";

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

    return passed;
}

bool testGetInt() 
{
    const std::string testIniFile = "test.ini";

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

    return passed;
}

bool testGetDouble(double absError) 
{   
    // Create an INI test file
    const std::string testIniFile = "test.ini";

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

    return passed;
}

bool testGetSectionData() 
{
    const std::string testIniFile = "test.ini";

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

    return passed;
}

bool testGetSections()
{
    const std::string testIniFile = "test.ini";

    // Create a test .ini file
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

    return passed;
}


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
