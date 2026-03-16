#!/bin/bash

# Build the test project using make
make -j$(nproc)

# Execute the tests
./test.out

# Clean up build artifacts using "make clean"
make clean
