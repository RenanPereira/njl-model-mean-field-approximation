#ifndef SU3NJL3DCUTOFFFILEPARSER_H
#define SU3NJL3DCUTOFFFILEPARSER_H

#include "ini_file_parser/IniFileParser.h"


class SU3NJL3DCutoffVacuumFileParser
{
public:
    const IniFileParser& parser;

    std::string sectionModelParameters = "SU3NJL3DCutoffModelParameters";
    
    std::string invalidFileMessage = "Error: Invalid configuration found in the " + parser.getFilename() + " file.";

public:
    SU3NJL3DCutoffVacuumFileParser(const IniFileParser& p) : parser(p) {}
    bool validateFileQuality() const;

};

#endif
