#ifndef INIFILEPARSER_H
#define INIFILEPARSER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>


class IniFileParser 
{
private:
    std::string filename;
    std::map<std::string, std::map<std::string, std::string>> sections;
    std::vector<std::string> sectionOrder; // Used to keep track of section order

public:
    IniFileParser(const std::string& );

public:
    std::string getFilename() const{ return filename; };

    // Method for getting a specific value from a given section
    std::string getValue(const std::string& , const std::string& ) const;

    // Method for getting all sections with a common prefix
    std::vector<std::map<std::string, std::string>> getSectionsData(const std::string& ) const;
    std::vector<std::pair<std::string, std::map<std::string, std::string>>> getSections(const std::string& ) const;

    // Method for getting an int variable
    int getInt(const std::string& , const std::string& ) const;
    int getInt(const std::string& , const std::string& , int ) const;
    int getInt(const std::map<std::string, std::string>& , const std::string& ) const;
    int getInt(const std::map<std::string, std::string>& , const std::string& , int ) const;

    // Method for getting a double variable
    double getDouble(const std::string& , const std::string& ) const;
    double getDouble(const std::string& , const std::string& , double ) const;
    double getDouble(const std::map<std::string, std::string>& , const std::string& ) const;
    double getDouble(const std::map<std::string, std::string>& , const std::string& , double ) const;

    static std::string trim(const std::string& );
};


#endif