#!/bin/bash

# This script runs tests for the modules which contain tests (akin to unitary tests and functional tests).

cd ../.. || exit

(cd src/ini_file_parser/ && ./execute_tests.sh)

(cd tests/integration_methods/ && ./execute_tests.sh)

(cd src/gsl_wrapper/ && ./execute_tests.sh)
