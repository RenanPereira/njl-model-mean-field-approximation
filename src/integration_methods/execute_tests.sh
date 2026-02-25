#!/bin/bash

# This script is designed to streamline the testing of the local module "integration_methods".
# Due to the structure of the project, it temporarily modifies the header file paths 
# in the source code to ensure compatibility with the 'make' build system. Once the 
# tests are executed, the script reverts all changes to restore the source files 
# to their original state.

# Modify "Integration1DimNewtonCotes.cpp" to update the include directive.
file="Integration1DimNewtonCotes.cpp"
search='#include "integration_methods/Integration1DimNewtonCotes.h"'
replace='#include "Integration1DimNewtonCotes.h"'
../../scripts/utils/modify_headers.sh "$file" "$search" "$replace"

# Build the test project using make
make

# Execute the tests
./tests.out

# Clean up build artifacts using "make clean"
make clean

# Revert the changes in "Integration1DimNewtonCotes.cpp" to its original state.
file="Integration1DimNewtonCotes.cpp"
../../scripts/utils/restore_headers.sh "$file"
