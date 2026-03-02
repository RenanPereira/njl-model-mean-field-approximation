#ifndef COMMAND_LINE_PROCESSOR_H
#define COMMAND_LINE_PROCESSOR_H

#include "ini_file_parser/IniFileParser.h"

int commandLineArgsProcessor(int , char* []);

void printQualityCheckFailedMessage(std::string );

void selectPathBasedOnFileDetails(const IniFileParser& );

#endif