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

    std::string getFilename() const{ return filename; };
    std::map<std::string, std::map<std::string, std::string>> getSections() const{ return sections; };
    std::vector<std::string> getSectionOrder() const{ return sectionOrder; };

    // Method for getting a specific value from a given section
    std::string getValue(const std::string& , const std::string& ) const;
    std::string getValue(const std::map<std::string, std::string>& , const std::string& ) const;

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

    bool isKeyPresent(const std::string& section, const std::string& key) const;
    
    bool validatePositiveInteger(const std::string& section, const std::string& key, 
                                 const std::string& invalidFileMessage, const std::string& conditionMessage) const;
    
    bool validatePositiveDouble(const std::string& section, const std::string& key, 
                                const std::string& invalidFileMessage, const std::string& conditionMessage) const;
    
    bool validateNonNegativeDouble(const std::string& section, const std::string& key, 
                                   const std::string& invalidFileMessage, const std::string& conditionMessage) const;
    
    bool validateRequiredSections(const std::vector<std::string>& ) const;
    
    bool validateRequiredKeys(const std::string& , const std::vector<std::string>& ) const;
    bool validateRequiredKeys(const std::map<std::string, std::string>& , const std::vector<std::string>& ) const;

    bool getBool(const std::string& , const std::string& , bool ) const;
    bool getBool(const std::string& , const std::string& ) const;
    bool validateBool(const std::string& , const std::string& , const std::string& , const std::string& ) const;
};


#endif