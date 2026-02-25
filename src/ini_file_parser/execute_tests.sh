#!/bin/bash

# This script is designed to streamline the testing of the local module "ini_file_parser".
# Due to the structure of the project, it temporarily modifies the header file paths 
# in the source code to ensure compatibility with the 'make' build system. Once the 
# tests are executed, the script reverts all changes to restore the source files 
# to their original state.

# Modify "IniFileParser.cpp" to update the include directive.
file="IniFileParser.cpp"
search='#include "ini_file_parser/IniFileParser.h"'
replace='#include "IniFileParser.h"'
../../scripts/utils/modify_headers.sh "$file" "$search" "$replace"

# Build the test project using make
make

# Execute the tests
./tests.out

# Clean up build artifacts using "make clean"
make clean

# Revert the changes in "IniFileParser.cpp" to its original state.
file="IniFileParser.cpp"
../../scripts/utils/restore_headers.sh "$file"
