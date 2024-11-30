#include "ini_file_parser/IniFileParser.h"


IniFileParser::IniFileParser(const std::string& filenameAux) 
{   
    // set the filename
    filename = filenameAux;

    std::ifstream file(filename);
    if (!file) 
	{
        std::cout << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string currentSection;

    // iterate over each line in the file, reading it into line, until the end of the file is reached
    while (std::getline(file, line)) 
	{   
        // remove any leading and trailing whitespace from the line
        line = trim(line);

        // check if the line is empty or if it is a comment
        if (line.empty() || line[0] == ';' || line[0] == '#') { continue; }


        // check if line corresponds to section, otherwise it should be a key=value pair, containing an "=" sign
        if (line[0] == '[') 
		{
            size_t closingBracket = line.find(']');

            // check if closing bracket is found in the string
            if (closingBracket == std::string::npos) //std::string::npos is returned when function cannot find the requested substring or character
			{
                std::cout << "Malformed section header: " << line << std::endl;
                abort();
            }

            // extract the text between [ and ] as the section name
            currentSection = trim(line.substr(1, closingBracket - 1));
            
            // checks if the section name is already in the sections map in order to prevent duplicates
            if (sections.find(currentSection) == sections.end()) 
            {
                sections[currentSection] = {};
                sectionOrder.push_back(currentSection); // Record section name in order
            }
            else
            {
                std::cout << "I detected that the section: " << currentSection << " is duplicated in the file: " << filename << "\n";
                std::cout << "Section names must be unique. Aborting.\n";
                abort();
            }
        } 
		else 
		{   
            // find the position of the "=" character which separates the key-value pair
            size_t pos = line.find('=');

            if (pos == std::string::npos) //std::string::npos is returned when function cannot find the requested substring or character
			{
                std::cout << "Malformed line (missing '=' sign): " << line << std::endl;
                abort();
            }

            // extract the key and value pairs (key is before "=", and value is after "=")
            std::string key = trim(line.substr(0, pos));
            std::string value = trim(line.substr(pos + 1));

            // add key-value pair if the section is not empty
            if (!currentSection.empty()) 
			{
                sections[currentSection][key] = value;
            }
        }
    }
}

std::string IniFileParser::getValue(const std::string& section, const std::string& key) const 
{
    auto secIt = sections.find(section);
    if (secIt != sections.end()) 
	{
        auto keyIt = secIt->second.find(key);
        if (keyIt != secIt->second.end()) 
    	{
            return keyIt->second;
        }
    }

    return "";
}

std::vector<std::map<std::string, std::string>> IniFileParser::getSectionsData(const std::string& sectionPrefix) const 
{
    std::vector<std::map<std::string, std::string>> results;
    for (int i = 0; i < int(sectionOrder.size()); ++i) // Use sectionOrder for correct ordering
    {
        const std::string& sectionName = sectionOrder[i];
        if (sectionName.find(sectionPrefix) == 0) 
        {
            results.push_back(sections.at(sectionName));
        }
    }

    return results;
}

std::vector<std::pair<std::string, std::map<std::string, std::string>>> IniFileParser::getSections(const std::string& sectionPrefix) const 
{
    std::vector<std::pair<std::string, std::map<std::string, std::string>>> results;

    for (int i = 0; i < int(sectionOrder.size()); ++i) // Use sectionOrder for correct ordering
    {   
        const std::string& sectionName = sectionOrder[i];
        if (sectionName.find(sectionPrefix) == 0) 
        {
            results.emplace_back(sectionName, sections.at(sectionName));
        }
    }

    return results;
}

int IniFileParser::getInt(const std::string& section, const std::string& key) const 
{
    std::string value = getValue(section, key);
    if (value.empty()) 
    {
        return 0;
    }
    else
    {
        return std::stoi(value);
    }
}

int IniFileParser::getInt(const std::string& section, const std::string& key, int defaultValue) const 
{
    std::string value = getValue(section, key);
    if (value.empty()) 
    {
        return defaultValue;
    }
    else
    {
        return std::stoi(value);
    }
}

int IniFileParser::getInt(const std::map<std::string, std::string>& section, const std::string& key) const 
{
    auto it = section.find(key);
    if ( it != section.end() )
    {
        return std::stoi(it->second);
    }
    else
    {
        return 0;
    }   
}

int IniFileParser::getInt(const std::map<std::string, std::string>& section, const std::string& key, int defaultValue) const 
{
    auto it = section.find(key);
    if ( it != section.end() )
    {
        return std::stoi(it->second);
    }
    else
    {
        return defaultValue;
    }   
}

double IniFileParser::getDouble(const std::string& section, const std::string& key) const 
{   
    std::string value = getValue(section, key);
    if (value.empty()) 
    {
        return 0.0;
    }
    else
    {
        return std::stod(value);
    }
}

double IniFileParser::getDouble(const std::string& section, const std::string& key, double defaultValue) const 
{
    std::string value = getValue(section, key);
    if (value.empty()) 
    {
        return defaultValue;
    }
    else
    {
        return std::stod(value);
    }
}

double IniFileParser::getDouble(const std::map<std::string, std::string>& section, const std::string& key) const 
{
    auto it = section.find(key);
    if ( it != section.end() )
    {
        return std::stod(it->second);
    }
    else
    {
        return 0.0;
    }   
}

double IniFileParser::getDouble(const std::map<std::string, std::string>& section, const std::string& key, double defaultValue) const 
{
    auto it = section.find(key);
    if ( it != section.end() )
    {
        return std::stod(it->second);
    }
    else
    {
        return defaultValue;
    }   
}

std::string IniFileParser::trim(const std::string& str) 
{
    const char* whitespace = " \t\n\r\f\v";
    size_t start = str.find_first_not_of(whitespace);
    size_t end = str.find_last_not_of(whitespace);
    if (start == std::string::npos || end == std::string::npos)
	{
		return "";
	}
	else
	{
		return str.substr(start, end - start + 1);
	}
}
