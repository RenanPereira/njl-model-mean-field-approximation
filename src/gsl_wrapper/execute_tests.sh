#!/bin/bash


# This script is designed to streamline the testing of the local module "gsl_wrapper".
# Due to the structure of the project, it temporarily modifies the header file paths 
# in the source code to ensure compatibility with the 'make' build system. Once the 
# tests are executed, the script reverts all changes to restore the source files 
# to their original state.


# Navigate to the "integration_methods" directory to modify files for testing and
# modify the header include path in Integration1DimNewtonCotes.cpp.
# After that, navigate back to the "gsl_wrapper" directory for further modifications.
cd .. && cd integration_methods
file="Integration1DimNewtonCotes.cpp"
search='#include "integration_methods/Integration1DimNewtonCotes.h"'
replace='#include "Integration1DimNewtonCotes.h"'
../../tools/module_tests/modify_headers.sh "$file" "$search" "$replace"
cd .. && cd gsl_wrapper

# Modify "Integration1DimGSL.cpp" to update the include directive.
file="Integration1DimGSL.cpp"
search='#include "gsl_wrapper/Integration1DimGSL.h"'
replace='#include "Integration1DimGSL.h"'
../../tools/module_tests/modify_headers.sh "$file" "$search" "$replace"

# Modify "Integration1DimGSL.h" to update the include directive.
file="Integration1DimGSL.h"
search='#include "integration_methods/Integration1DimNewtonCotes.h"'
replace='#include "Integration1DimNewtonCotes.h"'
../../tools/module_tests/modify_headers.sh "$file" "$search" "$replace"

# Build the test project using make
make

# Execute the tests
./tests.out

# Clean up build artifacts using "make clean"
make clean

# Revert the changes in "Integration1DimGSL.cpp" to its original state.
file="Integration1DimGSL.cpp"
../../tools/module_tests/restore_headers.sh "$file"

# Revert the changes in "Integration1DimGSL.h" to its original state.
file="Integration1DimGSL.h"
../../tools/module_tests/restore_headers.sh "$file"

# Navigate back to "integration_methods" and restore the original file and return to the "gsl_wrapper" directory
cd .. && cd integration_methods
file="Integration1DimNewtonCotes.cpp"
../../tools/module_tests/restore_headers.sh "$file"
cd .. && cd gsl_wrapper
